//! Native local-hardness Archard integration and original normal continuation.

use fs_exec::Cx;
use fs_tribo::{ArchardLaw, LocalHardnessArchardSample, WearAdvance};

use super::{
    HistoryHardnessContactError, HistoryHardnessControl, HistoryHardnessWearAdvance,
    SourceHistoryHardnessWear, checkpoint, invalid, within,
};

impl SourceHistoryHardnessWear {
    /// Remove material using each body's mapped accepted hardness and the actual
    /// unsigned cell slip. Beginning pressure and fixed k/H are held within this
    /// interval; refine slip when removal changes pressure appreciably. Neither
    /// removal nor reduced indentation advances the external bulk history.
    ///
    /// Only the complete native re-equilibration and both volume/profile accounts
    /// can be accepted. A late source, numerical, budget or cancellation refusal
    /// leaves the entire predecessor unchanged.
    pub fn advance_wear(
        &self,
        unsigned_slip_m: &[f64],
        control: HistoryHardnessControl,
        cx: &Cx<'_>,
    ) -> Result<HistoryHardnessWearAdvance, HistoryHardnessContactError> {
        self.admit_control(control, true, cx)?;
        let pressure = self.contact.pressure_pa();
        if unsigned_slip_m.len() != pressure.len()
            || unsigned_slip_m.iter().any(|s| !s.is_finite() || *s < 0.0)
        {
            return Err(invalid(
                "finite nonnegative actual unsigned slip at every original grid cell required",
            ));
        }
        let grid = self.model.config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
        let mut recession = self.recession_m.clone();
        let mut increment = [vec![0.0; pressure.len()], vec![0.0; pressure.len()]];
        let mut relative = vec![0.0; pressure.len()];
        let mut measure = 0.0;
        let mut measure_correction = 0.0;
        let mut geometric = [0.0; 2];
        let mut corrections = [0.0; 2];
        let mut roundoff = [0.0; 2];
        for (i, (&p, &slip)) in pressure.iter().zip(unsigned_slip_m).enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            let cell_measure = (p * area) * slip;
            if !cell_measure.is_finite() || (p > 0.0 && slip > 0.0 && cell_measure == 0.0) {
                return Err(invalid(
                    "local pressure-sliding measure lies outside arithmetic resolution",
                ));
            }
            add(&mut measure, &mut measure_correction, cell_measure);
            for side in 0..2 {
                let h = self.hardness[side]
                    .hardness_pa(self.sources.cell_to_material_point[i][side])
                    .ok_or_else(|| {
                        invalid("original material-point map no longer addresses its field")
                    })?;
                let depth = (self.sources.k[side] * (p / h)) * slip;
                let old = recession[side][i];
                let total = old + depth;
                let retained = total - old;
                let active = self.sources.k[side] > 0.0 && p > 0.0 && slip > 0.0;
                if !depth.is_finite() || !total.is_finite() || (active && retained == 0.0) {
                    return Err(invalid(
                        "local-hardness recession lies outside retained profile resolution",
                    ));
                }
                within(
                    "body cumulative recession [m]",
                    total,
                    self.sources.maximum_recession_m[side],
                )?;
                let volume = retained * area;
                if !volume.is_finite() || (retained > 0.0 && volume == 0.0) {
                    return Err(invalid(
                        "material removal lies outside volume arithmetic resolution",
                    ));
                }
                add(&mut geometric[side], &mut corrections[side], volume);
                roundoff[side] += 4.0 * f64::EPSILON * ((old * area).abs() + (total * area).abs());
                increment[side][i] = retained;
                recession[side][i] = total;
            }
            relative[i] = increment[0][i] + increment[1][i];
            if !relative[i].is_finite() {
                return Err(invalid("finite relative material removal required"));
            }
        }
        if !measure.is_finite() {
            return Err(invalid("finite total pressure-weighted sliding required"));
        }
        let mut volumes = self.volumes;
        let interface = self.sources.normal.elastic().interface().interface();
        let mut advance_volume =
            |side: usize| -> Result<WearAdvance, HistoryHardnessContactError> {
                let mut refused = None;
                let samples = pressure.iter().zip(unsigned_slip_m).enumerate().scan(
                    (),
                    |_, (i, (&p, &slip))| {
                        if i % 64 == 0 {
                            if let Err(error) = checkpoint(cx) {
                                refused = Some(error);
                                return None;
                            }
                        }
                        Some(LocalHardnessArchardSample {
                            wear_coefficient: self.sources.k[side],
                            pressure_weighted_sliding_j: (p * area) * slip,
                            hardness_pa: self.hardness[side]
                                .hardness_pa(self.sources.cell_to_material_point[i][side])
                                .expect("immutable mapping admitted at construction"),
                        })
                    },
                );
                let result = ArchardLaw::advance_local_hardness_pressure_weighted_sliding(
                    interface,
                    &mut volumes[side],
                    samples,
                );
                if let Some(error) = refused {
                    return Err(error);
                }
                Ok(result?)
            };
        let wear = [advance_volume(0)?, advance_volume(1)?];
        for side in 0..2 {
            let native = wear[side].increment_m3;
            let retained = volumes[side].volume_m3() - self.volumes[side].volume_m3();
            let arithmetic =
                roundoff[side] + 16.0 * f64::EPSILON * native.abs().max(geometric[side].abs());
            let inventory = arithmetic
                + 4.0
                    * f64::EPSILON
                    * (self.volumes[side].volume_m3().abs() + volumes[side].volume_m3().abs());
            if !arithmetic.is_finite()
                || !inventory.is_finite()
                || !geometric[side].is_finite()
                || (native > 0.0 && retained == 0.0)
                || (native - geometric[side]).abs() > arithmetic
                || (retained - geometric[side]).abs() > inventory
            {
                return Err(invalid(
                    "original per-body Archard inventory and actual recession disagree beyond roundoff",
                ));
            }
        }
        let (step, control) = self.recede(&relative, control, cx)?;
        self.admit(&self.model, &step.candidate, &recession, cx)?;
        let candidate = Self {
            sources: std::sync::Arc::clone(&self.sources),
            hardness: self.hardness.clone(),
            model: self.model.clone(),
            contact: step.candidate,
            volumes,
            recession_m: recession,
        };
        Ok(HistoryHardnessWearAdvance {
            candidate,
            wear,
            recession_increment_m: increment,
            pressure_weighted_sliding_j: measure,
            receipt: step.receipt,
            control,
        })
    }
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let corrected = value - *correction;
    let next = *sum + corrected;
    *correction = (next - *sum) - corrected;
    *sum = next;
}
