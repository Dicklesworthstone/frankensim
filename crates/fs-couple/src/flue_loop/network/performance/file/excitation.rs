//! One source-law selection for characteristic networks and spatial DG bores.
//!
//! The bore owner changes neither the windway's material law nor the jet's
//! admissible profile. These records attach the original owners in the same
//! order, retaining their incompatibility, gas and validity checks.

use super::{FlueInputError, Rows, bad, poll};
use crate::bernoulli_aperture::contraction::ContractionSpec;
use crate::flue_loop::auvray_jet::AuvraySpec;
use crate::flue_loop::forced_jet::{ForcedJetSpec, JetExitLaw, JetSpreadingLaw};
use crate::flue_loop::jet_swing_limit;
use crate::flue_loop::momentum_drive;
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

/// One `jet_transfer` row: the Nolle/Price transfer or Auvray's jet drive.
enum JetTransferInput {
    Nolle(NolleSpec),
    Auvray(AuvraySpec),
}

/// Optional source selections after `jet`, before `jet_limits`.
pub(crate) struct FlueJetInput {
    card_centreline: bool,
    transfer: Option<JetTransferInput>,
    response: Option<(ForcedJetSpec, JetExitLaw)>,
    momentum: Option<f64>,
    spreading: JetSpreadingLaw,
    swing_limit: bool,
    cut_width_station: bool,
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
            match row.fields.next() {
                Some("nolle-price-2015") => {
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
                    Some(JetTransferInput::Nolle(
                        NolleSpec::price_2015(channel, chamfer_depth_m)
                            .map_err(|what| FlueInputError::Physics(FlueError::Invalid { what }))?,
                    ))
                }
                Some("auvray-2014") => {
                    // `jet_transfer auvray-2014 BETA GAMMA SOURCE_TAG`: every field required, none defaulted.
                    let beta = row.scalar()?;
                    let gamma = row.scalar()?;
                    let source = row.fields.next().ok_or_else(|| {
                        bad(row.line, "jet_transfer auvray-2014 needs an explicit source tag after BETA GAMMA")
                    })?;
                    row.end()?;
                    Some(JetTransferInput::Auvray(
                        AuvraySpec::new(beta, gamma, source)
                            .map_err(|what| FlueInputError::Physics(FlueError::Invalid { what }))?,
                    ))
                }
                _ => {
                    return Err(bad(
                        row.line,
                        "jet_transfer requires nolle-price-2015 CHANNEL CHAMFER_DEPTH_M or auvray-2014 BETA GAMMA SOURCE_TAG",
                    ));
                }
            }
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
            let exit = match row.fields.next() {
                None => JetExitLaw::Kutta,
                Some("exit") => match row.fields.next() {
                    Some("kutta") => JetExitLaw::Kutta,
                    Some("measured") => JetExitLaw::Measured,
                    _ => {
                        return Err(bad(
                            row.line,
                            "jet_response exit must be kutta or measured",
                        ));
                    }
                },
                Some(_) => {
                    return Err(bad(
                        row.line,
                        "jet_response requires forced-bickley SPACING WINDOW [exit kutta|measured]",
                    ));
                }
            };
            row.end()?;
            Some((spec, exit))
        } else {
            None
        };
        // Optional `momentum_drive KAPPA` right after `jet_response`: KAPPA is required and only the two printed weights
        // (1, Thwaites and Fletcher 1983 Eq. (1); 2.6, their Sec. II fit) are accepted. It acts on the forced jet's cut flow.
        let momentum = if rows.next_is("momentum_drive") {
            let mut row = rows.row("momentum_drive")?;
            let line = row.line;
            if row.fields.clone().next().is_none() {
                return Err(bad(line, momentum_drive::KAPPA_REFUSED));
            }
            let kappa = row.scalar().map_err(|_| bad(line, momentum_drive::KAPPA_REFUSED))?;
            row.end()?;
            if response.is_none() {
                return Err(bad(line, "momentum_drive requires the jet_response row before it"));
            }
            if !momentum_drive::is_printed_weight(kappa) {
                return Err(bad(line, momentum_drive::KAPPA_REFUSED));
            }
            Some(kappa)
        } else {
            None
        };
        // Optional `jet_spreading LAW` after `jet_response` (and `momentum_drive`, when present): the
        // forced jet's base widens by the named law, `bickley` or `fletcher-1976`; no default law.
        let spreading = if rows.next_is("jet_spreading") {
            let mut row = rows.row("jet_spreading")?;
            let line = row.line;
            let law = match row.fields.next() {
                Some("bickley") => JetSpreadingLaw::Bickley,
                Some("fletcher-1976") => JetSpreadingLaw::Fletcher1976,
                _ => return Err(bad(line, "jet_spreading requires bickley or fletcher-1976")),
            };
            row.end()?;
            if response.is_none() {
                return Err(bad(line, "jet_spreading requires the jet_response row before it"));
            }
            if rows.next_is("momentum_drive") {
                return Err(bad(line + 1, "momentum_drive goes before jet_spreading"));
            }
            law
        } else {
            JetSpreadingLaw::None
        };
        // Optional `jet_swing_limit fletcher-1976` after `jet_spreading`: Fletcher's envelope on the
        // forced jet's swing. It takes the spread width at the lip, so it needs the law row.
        let swing_limit = rows.next_is("jet_swing_limit");
        if swing_limit {
            let mut row = rows.row("jet_swing_limit")?;
            let line = row.line;
            if row.fields.next() != Some(jet_swing_limit::NAME) {
                return Err(bad(line, "jet_swing_limit requires fletcher-1976"));
            }
            row.end()?;
            if spreading == JetSpreadingLaw::None {
                return Err(bad(line, "jet_swing_limit requires the jet_spreading row before it"));
            }
        }
        // Optional `jet_cut_width station`, last: the cut takes the forced jet's own width and the
        // printed swing.
        let cut_width_station = rows.next_is("jet_cut_width");
        if cut_width_station {
            let mut row = rows.row("jet_cut_width")?;
            let line = row.line;
            if row.fields.next() != Some("station") {
                return Err(bad(line, "jet_cut_width requires station"));
            }
            row.end()?;
            if response.is_none() {
                return Err(bad(line, "jet_cut_width requires the jet_response row before it"));
            }
            if rows.next_is("jet_swing_limit") {
                return Err(bad(line + 1, "jet_swing_limit goes before jet_cut_width"));
            }
        }
        if transfer.is_some() && response.is_some() {
            return Err(bad(
                rows.line,
                "jet_transfer and jet_response exclude each other",
            ));
        }
        // A second jet row in either place would otherwise surface as a missing `jet_limits`.
        if rows.next_is("jet_transfer") || rows.next_is("jet_response") {
            return Err(bad(
                rows.line + 1,
                "one jet_transfer or jet_response row at most: the jet transfers exclude each other",
            ));
        }
        Ok(Self {
            card_centreline,
            transfer,
            response,
            momentum,
            spreading,
            swing_limit,
            cut_width_station,
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
            Some(JetTransferInput::Nolle(spec)) => {
                poll(gate)?;
                let voice = voice
                    .with_nolle_jet(spec)
                    .map_err(FlueInputError::Physics)?;
                poll(gate)?;
                voice
            }
            Some(JetTransferInput::Auvray(spec)) => {
                poll(gate)?;
                let voice = voice
                    .with_auvray_jet(spec)
                    .map_err(FlueInputError::Physics)?;
                poll(gate)?;
                voice
            }
            None => voice,
        };
        let voice = match self.response {
            Some((spec, exit)) => {
                poll(gate)?;
                let voice = voice
                    .with_forced_jet_spreading(spec, exit, self.spreading, held_pressure_pa)
                    .map_err(FlueInputError::Physics)?;
                poll(gate)?;
                voice
            }
            None => voice,
        };
        let voice = match self.momentum {
            Some(kappa) => voice
                .with_momentum_drive(kappa)
                .map_err(FlueInputError::Physics)?,
            None => voice,
        };
        let voice = if self.swing_limit {
            voice.with_jet_swing_limit().map_err(FlueInputError::Physics)?
        } else {
            voice
        };
        if self.cut_width_station {
            voice.with_jet_cut_width_station().map_err(FlueInputError::Physics)
        } else {
            Ok(voice)
        }
    }
}
