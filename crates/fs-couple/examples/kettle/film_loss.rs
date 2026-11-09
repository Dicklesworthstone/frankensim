//! Opt-in head film loss (`--head-film-loss PACK`): a canonical material card's
//! tensile loss-tangent curve, put on each head mode's BENDING share only.
//! The film carries its installed tension as prestress, and prestress is not
//! lossy strain, so the curve never multiplies the tension's stiffness: on
//! this simply supported head (w = 0, free slope) the bending share is the
//! mode's energy in the same DKT pencil assembled with no prestress. A mode's
//! damping in the head's `damping_per_s` convention is then
//! `tan_delta(f) * share * omega`, with f and omega its in-vacuo values,
//! a Kelvin-Voigt bending law that matches the curve at each mode's own
//! frequency. The card's curve was measured on another film than this head
//! (its own receipt says which); the line printed per mode says so.
use super::Error;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_plate::{AssemblyOptions, EdgeSupport, ModePair, assemble, shell::head::TensionedDisk};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use std::{f64::consts::TAU, io::Read};

pub const FLAG: &str = "--head-film-loss";
/// The card property read, on its typed cyclic `frequency` axis.
pub const PROPERTY: &str = "loss_tangent_tensile";
pub const TEMPERATURE_AXIS: &str = "T";
pub const FREQUENCY_AXIS: &str = "frequency";
/// The film's temperature: the card curve's only valid temperature, and the
/// kit's standard gas temperature.
pub const TEMPERATURE_K: f64 = 293.15;
const MAX_PACK_BYTES: u64 = 16 << 20;
pub const VISCOUS_RATIO_REFUSAL: &str =
    "--head-film-loss replaces the card's uniform viscous_ratio, which damps the prestress too; a card with a nonzero viscous_ratio refuses the flag";

pub struct FilmLoss {
    pack: NormalizedMaterialCardPack,
    path: String,
}
impl std::fmt::Debug for FilmLoss {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FilmLoss").field("path", &self.path).finish_non_exhaustive()
    }
}

/// Removes `--head-film-loss PACK` from the arguments and reads the pack.
/// # Errors
/// A missing or repeated path, an unreadable or oversized pack, or a pack
/// that does not decode.
pub fn option(args: &mut Vec<String>) -> Result<Option<FilmLoss>, Error> {
    let Some(i) = args.iter().position(|a| a == FLAG) else { return Ok(None) };
    if i + 1 >= args.len() { return Err(format!("{FLAG} needs a material card pack path").into()); }
    let path = args.remove(i + 1);
    args.remove(i);
    if args.iter().any(|a| a == FLAG) { return Err(format!("{FLAG} is given more than once").into()); }
    FilmLoss::load(&path).map(Some)
}

impl FilmLoss {
    /// # Errors
    /// An unreadable or oversized pack, or one that does not decode.
    pub fn load(path: &str) -> Result<Self, Error> {
        let mut bytes = Vec::new();
        std::fs::File::open(path)?.take(MAX_PACK_BYTES + 1).read_to_end(&mut bytes)?;
        if bytes.len() as u64 > MAX_PACK_BYTES { return Err(format!("{FLAG}: the pack exceeds 16 MiB").into()); }
        let pack = NormalizedMaterialCardPack::from_bytes(&bytes).map_err(|e| format!("{FLAG}: {e}"))?;
        Ok(Self { pack, path: path.into() })
    }

    /// The card's loss tangent at `frequency_hz` and TEMPERATURE_K, through the
    /// canonical resolver (its validity and its interpolation between knots).
    /// # Errors
    /// A frequency or temperature outside the claim's validity, or a card
    /// without exactly one such claim.
    pub fn tan_delta(&self, frequency_hz: f64) -> Result<f64, Error> {
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(format!("{FLAG}: frequency {frequency_hz} Hz").into());
        }
        let point = QueryPoint::new()
            .with_quantity(TEMPERATURE_AXIS,
                QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static)),
                TEMPERATURE_K)
            .and_then(|p| p.with_quantity(FREQUENCY_AXIS,
                QuantitySpec::semantic(SemanticType::new(QuantityKind::Frequency(FrequencyConvention::Cyclic), ValueForm::Static)),
                frequency_hz))
            .map_err(|e| e.to_string())?;
        let requirement = ScalarPropertyRequirement::try_new(PROPERTY, Dims::NONE, ScalarAdmissibility::StrictlyPositive)
            .map_err(|e| e.to_string())?;
        let value = resolve_material_state_point(self.pack.card(), &point, &[requirement],
            MaterialPropertySelection::SingleClaimOnly)
            .map_err(|e| format!("{FLAG}: the card does not resolve {PROPERTY} at {frequency_hz} Hz: {e}"))?
            .property(PROPERTY)
            .map(|p| p.value_si())
            .ok_or_else(|| format!("{FLAG}: {PROPERTY} did not resolve"))?;
        if !value.is_finite() || value <= 0.0 { return Err(format!("{FLAG}: {PROPERTY} {value}").into()); }
        Ok(value)
    }

    /// One damping rate per head mode [1/s, the head's `damping_per_s`
    /// convention], with one receipt line per mode.
    /// # Errors
    /// See [`bending_shares`] and [`Self::tan_delta`].
    pub fn damping(&self, film: &TensionedDisk, modes: &[ModePair], rim_nodes: &[usize])
        -> Result<(Vec<f64>, Vec<String>), Error> {
        let shares = bending_shares(film, modes, rim_nodes)?;
        let mut damping = Vec::with_capacity(modes.len());
        let mut lines = Vec::with_capacity(modes.len());
        for (n, (mode, share)) in modes.iter().zip(shares).enumerate() {
            let omega = mode.lambda.sqrt();
            let frequency_hz = omega / TAU;
            let tan_delta = self.tan_delta(frequency_hz)?;
            let rate = tan_delta * share * omega;
            let sigma = rate / 2.0;
            lines.push(format!(
                "head film loss (opt-in {FLAG} {}): mode {n} in vacuo {frequency_hz:.3} Hz, bending share {share:.6e}, \
                 tan delta {tan_delta:.6} at {TEMPERATURE_K} K, sigma {sigma:.6e} 1/s ({:.6} dB/s); on bending only \
                 (prestress carries no loss); the card's curve is measured on another film than this head",
                self.path, 20.0 / std::f64::consts::LN_10 * sigma));
            damping.push(rate);
        }
        Ok((damping, lines))
    }
}

/// Each mode's bending share `phi^T K_b phi / lambda`: K_b is the SAME mesh
/// and section assembled with no prestress, simply supported on the same rim,
/// so `phi` (mass-normalized) applies unchanged.
/// # Errors
/// An assembly that does not keep the film's reduced coordinates, or a share
/// outside [0, 1).
pub fn bending_shares(film: &TensionedDisk, modes: &[ModePair], rim_nodes: &[usize]) -> Result<Vec<f64>, Error> {
    let bending = assemble(&film.mesh, &film.section, rim_nodes, &[],
        &AssemblyOptions { pretension: 0.0, support: EdgeSupport::SimplySupported })?;
    if bending.free != film.model.free || bending.dof_map != film.model.dof_map {
        return Err("bending-only assembly does not keep the film's reduced coordinates".into());
    }
    let mut product = vec![0.0; bending.free];
    modes.iter().map(|mode| {
        if mode.phi.len() != bending.free || !mode.lambda.is_finite() || mode.lambda <= 0.0 {
            return Err::<f64, Error>("head mode does not match the film's pencil".into());
        }
        bending.k.spmv(&mode.phi, &mut product);
        let share = mode.phi.iter().zip(&product).map(|(a, b)| a * b).sum::<f64>() / mode.lambda;
        if !share.is_finite() || !(0.0..1.0).contains(&share) {
            return Err(format!("bending share {share} is outside [0, 1)").into());
        }
        Ok(share)
    }).collect()
}
