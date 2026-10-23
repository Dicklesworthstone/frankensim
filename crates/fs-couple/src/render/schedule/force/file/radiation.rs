//! Observer transfers derived from declared modal geometry and the surrounding gas.
//!
//! An optional block after the last force event replaces every authored
//! narrow-band transfer of the file:
//! ```text
//! radiation compact-baffled
//! ambient GAS_RECORD
//! observer DISTANCE_M
//! areas MODE_COUNT A_0 ... A_{MODE_COUNT-1}
//! ... exactly one areas row per voice, in voice order ...
//! ```
//! `A_k` is the signed modal volume `integral(phi_k dA)` [m^2/sqrt(kg)] of the
//! baffled radiating face in the file's mass-normalized basis. The gas record
//! is the shared `read_gas_record` declaration. Each mode then observes
//! `H_k = -i omega_k rho A_k exp(i omega_k r / c) / (2 pi r)` under exp(-i omega t).
//! With the runtime's narrow-band realization `p = Re(H) qdot + Im(H) omega q`,
//! that is exactly `rho A_k qddot(t - r/c) / (2 pi r)` for a steady sinusoid at
//! the mode's natural frequency: density scales it, and the sound speed sets
//! its phase lag. It is still the narrow-band observation. Transients are not
//! retarded, and no reaction load, absorption or directivity is added. Mode
//! rows must carry zero H fields, so no authored gain survives. A free
//! coordinate has no pressure transfer and must declare a zero area.
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources, read_gas_record};
use fs_math::{c64::C64, det};
use super::{MAX_MODES, MAX_VOICES, ModalPerformanceError, Reader, input};

/// Exact first record of the block.
pub const COMPACT_RADIATION_RECORD: &str = "radiation compact-baffled";

/// Resolved observer medium and distance of a radiating modal file.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactRadiationInfo {
    /// Surrounding gas the transfers were computed in.
    pub gas: GasState,
    /// Observer distance from the baffled face [m].
    pub listener_m: f64,
}

pub(super) struct CompactRadiation {
    pub(super) gas_record: GasRecord,
    pub(super) listener_m: f64,
    areas: Vec<Vec<f64>>,
}

impl CompactRadiation {
    /// Byte offset of the newline that precedes the block, if the file has one.
    pub(super) fn offset(text: &str) -> Option<usize> { text.find("\nradiation ") }

    /// Read the block that starts after `text[..offset]`; line numbers continue
    /// from the mechanical records. The gas never defaults to air.
    pub(super) fn read(
        text: &str,
        offset: usize,
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, ModalPerformanceError> {
        let mut reader = Reader { lines: text[offset + 1..].lines(), line: text[..offset].lines().count() };
        let mut row = reader.row("radiation")?;
        if row.word()? != "compact-baffled" {
            return Err(input(row.line, "radiation model must be compact-baffled"));
        }
        row.finish()?;
        let mut row = reader.row("ambient")?;
        let line = row.line;
        let gas_record = read_gas_record(&mut row.fields, None, sources).map_err(|error| match error {
            GasInputError::Input(what) => input(line, what),
            GasInputError::Material(source) => ModalPerformanceError::Gas { line, source },
        })?;
        row.finish()?;
        let listener_m: f64 = reader.one("observer")?;
        if !(listener_m.is_finite() && listener_m > 0.0) {
            return Err(input(reader.line, "observer distance must be finite and positive"));
        }
        let mut areas = Vec::new();
        let mut total = 0;
        while reader.lines.clone().next().is_some() {
            if areas.len() == MAX_VOICES {
                return Err(input(reader.line + 1, "more areas rows than the voice budget"));
            }
            let mut row = reader.row("areas")?;
            let count = row.count(MAX_MODES - total)?;
            total += count;
            let mut values = Vec::with_capacity(count);
            for _ in 0..count { values.push(row.scalar()?); }
            row.finish()?;
            areas.push(values);
        }
        Ok(Self { gas_record, listener_m, areas })
    }

    pub(super) fn info(&self) -> CompactRadiationInfo {
        CompactRadiationInfo { gas: *self.gas_record.state(), listener_m: self.listener_m }
    }

    /// The transfer of mode `mode` of voice `voice` at its natural frequency.
    /// Free coordinates (omega = 0) require a zero area and observe nothing.
    pub(super) fn transfer(&self, voice: usize, mode: usize, omega: f64, line: usize)
        -> Result<C64, ModalPerformanceError>
    {
        let area = *self.areas.get(voice).and_then(|row| row.get(mode))
            .ok_or_else(|| input(line, "radiation areas row is missing this mode"))?;
        if omega == 0.0 {
            return if area == 0.0 { Ok(C64::ZERO) }
                else { Err(input(line, "a free coordinate has no pressure transfer; declare zero area")) };
        }
        Ok(compact_baffled_transfer(omega, area, self.gas_record.state(), self.listener_m))
    }

    /// Every voice has exactly one row, of exactly its mode count.
    pub(super) fn check_counts(&self, modes_per_voice: &[usize], line: usize)
        -> Result<(), ModalPerformanceError>
    {
        if self.areas.len() != modes_per_voice.len()
            || self.areas.iter().zip(modes_per_voice).any(|(row, &modes)| row.len() != modes)
        {
            return Err(input(line, "radiation needs one areas row per voice, each of its mode count"));
        }
        Ok(())
    }
}

/// `-i omega rho A exp(i omega r / c) / (2 pi r)`: a baffled compact monopole of
/// signed volume `A` per unit modal displacement, heard at distance `r`.
#[must_use]
pub fn compact_baffled_transfer(omega: f64, area: f64, gas: &GasState, listener_m: f64) -> C64 {
    let scale = omega * gas.density * area / (core::f64::consts::TAU * listener_m);
    let phase = omega * listener_m / gas.sound_speed;
    // -i (cos + i sin) = sin - i cos.
    C64::new(scale * det::sin(phase), -scale * det::cos(phase))
}
