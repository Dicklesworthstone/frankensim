//! Reciprocal diffuse sheets in the ORIGINAL participating enclosure.
//!
//! Each pair is one isothermal sheet with two equal-area faces. Its sourced
//! absorptance/emissivity and transmission belong to the conduction owner.
//! The original reciprocal chords still determine irradiation on each side:
//! `J_i = a_i B + r_i H_i + t H_partner`. Transmitted power is not heat;
//! only `A_i a_i (H_i - B)` enters that face's thermal account. Apply both
//! face powers to the same sheet thermal owner exactly once.
//!
//! This diffuse, equal-external-index boundary has no specular, refractive,
//! angle-resolved or finite-thickness temperature-gradient interpretation.

use fs_conduction::radiation::sheet::ReciprocalSheetOptics;

use super::{
    Layer, LayerAbsorption, MAX_RADIATION_SURFACES, ParticipatingEnclosure, RadiationError,
    RadiationReport, RadiationSolve, RadiationWorkspace, SurfaceEmissivity, Temperature, TwoZone,
    checkpoint, finite,
};

/// Two boundary faces of one sourced, reciprocal, isothermal diffuse sheet.
#[derive(Clone, Debug, PartialEq)]
pub struct SheetPair {
    /// Distinct indices in the original view-factor axis, in source-side order.
    pub surfaces: [usize; 2],
    /// Original absorptance, transmission, validity and source for both sides.
    pub optics: ReciprocalSheetOptics,
}

/// Admit the pairing/source axes without inferring geometry or temperatures.
pub(super) fn admit_axis(
    n: usize,
    emissivities: Option<&[f64]>,
    sheets: &[SheetPair],
) -> Result<[bool; MAX_RADIATION_SURFACES], RadiationError> {
    if n > MAX_RADIATION_SURFACES || sheets.len() > n / 2 {
        return Err(RadiationError::Invalid("invalid sheet surface axis"));
    }
    let mut paired = [false; MAX_RADIATION_SURFACES];
    for sheet in sheets {
        let [i, j] = sheet.surfaces;
        if i >= n || j >= n || i == j || paired[i] || paired[j] {
            return Err(RadiationError::Invalid(
                "sheet faces must be distinct, in range and belong to one pair",
            ));
        }
        if let Some(values) = emissivities {
            let a = sheet.optics.absorptance();
            if values.len() != n || values[i] != a[0] || values[j] != a[1] {
                return Err(RadiationError::Invalid(
                    "sheet absorptance must equal the original band emissivity on both faces",
                ));
            }
        }
        paired[i] = true;
        paired[j] = true;
    }
    Ok(paired)
}

impl ParticipatingEnclosure {
    /// Evaluate gray sheets and opaque patches with the same gas/path operator.
    ///
    /// `emissivities` contains ONLY unpaired opaque faces, in increasing
    /// view-factor order. Paired faces use their original sheet source, allowing
    /// exactly transparent/reflecting limits without inventing opaque cards.
    /// Equal paired area and temperature are required on every evaluation.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_with_sheets(
        &self,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        sheets: &[SheetPair],
        gas_temperatures: [Temperature; 2],
        absorption: &[LayerAbsorption; 2],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RadiationReport, RadiationError> {
        self.evaluate_gray(
            surface_temperatures,
            emissivities,
            sheets,
            gas_temperatures,
            absorption,
            controls,
            work,
            cancelled,
        )
    }

    /// Use original gas-layer temperatures and reject paths through absent gas.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_room_with_sheets(
        &self,
        room: &TwoZone,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        sheets: &[SheetPair],
        absorption: &[LayerAbsorption; 2],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiationReport, RadiationError> {
        for path in &self.paths {
            checkpoint(&mut cancelled)?;
            for segment in &path.segments {
                if segment.length.value() > 0.0 && room.layer(segment.layer).volume().value() == 0.0
                {
                    return Err(RadiationError::Invalid(
                        "optical path crosses an absent gas layer",
                    ));
                }
            }
        }
        self.evaluate_with_sheets(
            surface_temperatures,
            emissivities,
            sheets,
            [
                room.layer(Layer::Lower).temperature(),
                room.layer(Layer::Upper).temperature(),
            ],
            absorption,
            controls,
            work,
            cancelled,
        )
    }

    pub(super) fn admit_sheet_state(
        &self,
        sheets: &[SheetPair],
        surface_temperatures: &[Temperature],
        emissivities: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[bool; MAX_RADIATION_SURFACES], RadiationError> {
        let paired = admit_axis(self.geometry.len(), emissivities, sheets)?;
        for sheet in sheets {
            checkpoint(cancelled)?;
            let [i, j] = sheet.surfaces;
            let temperature = surface_temperatures[i].value();
            if self.geometry.areas_m2()[i] != self.geometry.areas_m2()[j]
                || temperature != surface_temperatures[j].value()
            {
                return Err(RadiationError::Invalid(
                    "paired sheet faces need equal physical areas and one isothermal temperature",
                ));
            }
            sheet
                .optics
                .validate_temperature(temperature)
                .map_err(|_| {
                    RadiationError::Invalid("sheet temperature is outside original optical support")
                })?;
        }
        Ok(paired)
    }

    pub(super) fn gray_sheet_emissivities(
        &self,
        temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        sheets: &[SheetPair],
        values: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RadiationError> {
        let paired = self.admit_sheet_state(sheets, temperatures, None, cancelled)?;
        if emissivities.len() != self.geometry.len() - 2 * sheets.len() {
            return Err(RadiationError::Invalid(
                "supply one original gray emissivity for each unpaired opaque face",
            ));
        }
        let mut sources = emissivities.iter();
        for (i, value) in values.iter_mut().enumerate() {
            if !paired[i] {
                let source = sources
                    .next()
                    .ok_or(RadiationError::Invalid("missing opaque source"))?;
                if source
                    .temperature_validity_k()
                    .is_some_and(|(lo, hi)| !(lo..=hi).contains(&temperatures[i].value()))
                {
                    return Err(RadiationError::Invalid(
                        "surface temperature is outside original emissivity validity",
                    ));
                }
                *value = source.value();
            }
        }
        for sheet in sheets {
            for (i, a) in sheet.surfaces.into_iter().zip(sheet.optics.absorptance()) {
                values[i] = a;
            }
        }
        Ok(())
    }

    // The opaque fast path retains the original matrix and arithmetic.
    pub(super) fn prepare_sheet_optics(
        &self,
        emissivities: &[f64],
        sheets: &[SheetPair],
        absorption: [f64; 2],
        work: &mut RadiationWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RadiationError> {
        if sheets.is_empty() {
            return self.prepare_optics(emissivities, absorption, work, cancelled);
        }
        self.prepare_exchange(absorption, work, cancelled)?;
        let n = self.geometry.len();
        let m = n + 2;
        let (reflection, transmission, partner) = coefficients(emissivities, sheets);
        work.matrix.fill(0.0);
        let mut anchored = [false; MAX_RADIATION_SURFACES];
        for i in 0..n {
            checkpoint(cancelled)?;
            let t = transmission[i];
            let p = partner[i];
            work.matrix[i * n + i] = emissivities[i] + t;
            work.matrix[i * n + p] -= t;
            anchored[i] = emissivities[i] > 0.0;
            // H_k = J_k + P_k/A_k. Keep the original implicit self-return
            // instead of assuming that the off-diagonal G row sums to A_k.
            for (k, fraction) in [(i, reflection[i]), (p, t)] {
                if fraction == 0.0 {
                    continue;
                }
                let scale = fraction / self.geometry.areas_m2()[k];
                for j in 0..m {
                    let c = finite(scale * work.exchange[k * m + j])?;
                    work.matrix[i * n + k] = finite(work.matrix[i * n + k] + c)?;
                    if j < n {
                        work.matrix[i * n + j] = finite(work.matrix[i * n + j] - c)?;
                    } else {
                        anchored[i] |= c > 0.0;
                    }
                }
            }
        }
        // Reverse reachability in the ACTUAL directed radiosity dependency
        // graph. Transmission can move an irradiation dependency across a
        // pair; an undirected chord-component check would miss that change.
        for _ in 0..n {
            checkpoint(cancelled)?;
            let mut changed = false;
            for i in 0..n {
                if !anchored[i] && (0..n).any(|j| anchored[j] && work.matrix[i * n + j] < 0.0) {
                    anchored[i] = true;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        if anchored[..n].iter().any(|&a| !a) {
            return Err(RadiationError::Invalid(
                "closed lossless sheet/surface component has indeterminate radiosity",
            ));
        }
        work.lu.factor(&work.matrix)?;
        checkpoint(cancelled)
    }

    pub(super) fn solve_sheet_emission(
        &self,
        emissivities: &[f64],
        sheets: &[SheetPair],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), RadiationError> {
        if sheets.is_empty() {
            return self.solve_emission(emissivities, controls, work, cancelled);
        }
        let n = self.geometry.len();
        let m = n + 2;
        let (reflection, transmission, partner) = coefficients(emissivities, sheets);
        for i in 0..n {
            checkpoint(cancelled)?;
            let mut rhs = finite(emissivities[i] * work.emission[i])?;
            for (k, fraction) in [(i, reflection[i]), (partner[i], transmission[i])] {
                if fraction == 0.0 {
                    continue;
                }
                let scale = fraction / self.geometry.areas_m2()[k];
                for g in n..m {
                    let c = finite(scale * work.exchange[k * m + g])?;
                    rhs = finite(c.mul_add(work.emission[g], rhs))?;
                }
            }
            work.rhs[i] = rhs;
        }
        work.lu
            .solve_factored_into(&work.rhs, &mut work.radiosity)?;
        checkpoint(cancelled)?;
        self.exchange_powers(work, cancelled)?;
        let mut irradiation = [0.0; MAX_RADIATION_SURFACES];
        for i in 0..n {
            irradiation[i] =
                finite(work.radiosity[i] + work.powers[i] / self.geometry.areas_m2()[i])?;
        }
        let mut maximum = 0.0_f64;
        for i in 0..n {
            let residual = finite(
                emissivities[i] * (work.radiosity[i] - work.emission[i])
                    - reflection[i] * work.powers[i] / self.geometry.areas_m2()[i]
                    + transmission[i] * (work.radiosity[i] - irradiation[partner[i]]),
            )?;
            maximum = maximum.max(residual.abs());
        }
        // Raw P includes through-transmission. Replace only sheet face heat;
        // retain gas/opaque exchanges and never force a balancing correction.
        for sheet in sheets {
            for (i, a) in sheet.surfaces.into_iter().zip(sheet.optics.absorptance()) {
                work.powers[i] = if a == 0.0 {
                    0.0
                } else {
                    finite(
                        (self.geometry.areas_m2()[i] * a)
                            * (work.radiosity[i] - work.emission[i]
                                + work.powers[i] / self.geometry.areas_m2()[i]),
                    )?
                };
            }
        }
        let balance = finite(work.powers.iter().copied().sum())?;
        if maximum > controls.surface_residual_w_m2 || balance.abs() > controls.energy_residual_w {
            return Err(RadiationError::Residual {
                surface_w_m2: maximum,
                energy_w: balance,
            });
        }
        checkpoint(cancelled)?;
        Ok((maximum, balance))
    }
}

fn coefficients(
    emissivities: &[f64],
    sheets: &[SheetPair],
) -> (
    [f64; MAX_RADIATION_SURFACES],
    [f64; MAX_RADIATION_SURFACES],
    [usize; MAX_RADIATION_SURFACES],
) {
    let mut reflection = [0.0; MAX_RADIATION_SURFACES];
    let mut transmission = [0.0; MAX_RADIATION_SURFACES];
    let mut partner = [0; MAX_RADIATION_SURFACES];
    for (i, &a) in emissivities.iter().enumerate() {
        reflection[i] = 1.0 - a;
        partner[i] = i;
    }
    for sheet in sheets {
        let [i, j] = sheet.surfaces;
        let [ri, rj] = sheet.optics.reflectance();
        reflection[i] = ri;
        reflection[j] = rj;
        transmission[i] = sheet.optics.transmittance();
        transmission[j] = sheet.optics.transmittance();
        partner[i] = j;
        partner[j] = i;
    }
    (reflection, transmission, partner)
}

#[cfg(test)]
mod tests;
