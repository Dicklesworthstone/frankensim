//! The plate-valve grammar consumes the shared material/geometry DG source.
//! Every bore and wall state stays in the original propagation owner.

use super::*;
use crate::dg_port::input::{DgSourceError, DgSourceInput};

pub(in super::super) type DgInput = DgSourceInput;

impl From<DgSourceError> for PlateValveInputError {
    fn from(error: DgSourceError) -> Self {
        match error {
            DgSourceError::Input { line, what } => Self::Input { line, what },
            DgSourceError::Gas { line, source } => Self::Gas { line, source },
            DgSourceError::Material { line, source } => Self::Material { line, source },
            DgSourceError::Physics(source) => Self::Physics(source),
            DgSourceError::Radiation { line, source } => Self::Physics(
                crate::acoustic_realize::AcousticRealizeError::Nonlinear(
                    format!("DG radiation line {line}: {source}"),
                ),
            ),
        }
    }
}
