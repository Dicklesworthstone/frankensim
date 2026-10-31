//! Phase-sensitive field recovery and spatial heat from the original sources.

use super::*;
use crate::planar::gradients;

impl HarmonicWindingFrequency<'_, '_, '_> {
    /// Superpose peak AC winding currents with arbitrary relative phases.
    /// Copper and resolved eddy Joule powers are reconstructed independently
    /// from their actual spatial sources; neither is inferred from terminal Z.
    pub fn response(
        &self,
        currents_a: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicWindingResponse> {
        poll(cancelled)?;
        let model = self.prepared;
        let count = model.basis.winding_count();
        if currents_a.len() != count
            || currents_a
                .iter()
                .any(|a| !a.re.is_finite() || !a.im.is_finite())
        {
            return Err(HarmonicWindingError::Invalid(
                "one finite peak phasor per original winding required",
            ));
        }
        let mesh = model.basis.operating_field().problem().mesh();
        let n = mesh.vertices().len();
        let mut potential = vec![C64::ZERO; n];
        let mut voltage = vec![C64::ZERO; count];
        for (j, (&current, response)) in currents_a.iter().zip(&self.unit_responses).enumerate() {
            for (value, &a) in potential
                .iter_mut()
                .zip(&response.potential_per_current_h_m)
            {
                poll(cancelled)?;
                *value = phasor(*value + a * current)?;
            }
            for (i, value) in voltage.iter_mut().enumerate() {
                *value = phasor(*value + self.impedance_ohm[i * count + j] * current)?;
            }
        }
        let mut flux_density = Vec::with_capacity(mesh.faces().len());
        for face in mesh.faces() {
            poll(cancelled)?;
            let shape = gradients(face.map(|v| mesh.vertices()[v as usize]))?;
            let a = face.map(|v| potential[v as usize]);
            let gradient: [C64; 2] = [0, 1].map(|axis| {
                (a[1] - a[0]).scale(shape[1][axis]) + (a[2] - a[0]).scale(shape[2][axis])
            });
            flux_density.push([phasor(gradient[1])?, phasor(-gradient[0])?]);
        }
        let real: Vec<_> = potential.iter().map(|a| a.re).collect();
        let imaginary: Vec<_> = potential.iter().map(|a| a.im).collect();
        let mut applied = vec![0.0; n];
        let mut mean_magnetic_energy_j = 0.0;
        for values in [&real, &imaginary] {
            poll(cancelled)?;
            model.tangent.spmv(values, &mut applied);
            for (&a, &ka) in values.iter().zip(&applied) {
                poll(cancelled)?;
                mean_magnetic_energy_j = finite((0.25 * a).mul_add(ka, mean_magnetic_energy_j))?;
            }
        }
        for (&current, &end) in currents_a.iter().zip(model.basis.end_inductance_h()) {
            mean_magnetic_energy_j =
                finite(mean_magnetic_energy_j + 0.25 * end * current.norm_sq())?;
        }
        if mean_magnetic_energy_j < 0.0 {
            return Err(HarmonicWindingError::Balance(
                "negative magnetic perturbation storage",
            ));
        }
        let conductors = solid::recover(
            model,
            &real,
            &imaginary,
            currents_a,
            self.omega_rad_s,
            cancelled,
        )?;
        let mut winding_cells = Vec::with_capacity(model.winding_conductors.len());
        let mut mean_winding_power_w = 0.0;
        for (source, &winding) in model
            .winding_conductors
            .iter()
            .zip(&model.conductor_windings)
        {
            poll(cancelled)?;
            let profile = model.basis.profiles_per_m2()[winding][source.triangle];
            let current_density_a_m2 = phasor(currents_a[winding].scale(profile))?;
            let electric_field_v_m =
                phasor(current_density_a_m2.scale(source.conductor.resistivity_ohm_m()))?;
            let volume_m3 =
                mesh.face_measure(source.triangle)
                    .ok_or(HarmonicWindingError::Invalid(
                        "original winding volume missing",
                    ))?;
            let mean_power_w = finite(
                0.5 * volume_m3
                    * (electric_field_v_m.re * current_density_a_m2.re
                        + electric_field_v_m.im * current_density_a_m2.im),
            )?;
            if mean_power_w < 0.0 || current_density_a_m2 != C64::ZERO && mean_power_w == 0.0 {
                return Err(HarmonicWindingError::Balance(
                    "represented winding heat is not positive",
                ));
            }
            mean_winding_power_w = finite(mean_winding_power_w + mean_power_w)?;
            winding_cells.push(HarmonicWindingCell {
                source: source.clone(),
                winding,
                volume_m3,
                current_density_a_m2,
                electric_field_v_m,
                mean_power_w,
            });
        }
        let mut supplied = C64::ZERO;
        for (&current, &v) in currents_a.iter().zip(&voltage) {
            supplied = phasor(supplied + (current.conj() * v).scale(0.5))?;
        }
        let power_residual_w = finite(
            conductors.mean_eddy_power_w + conductors.mean_solid_power_w + mean_winding_power_w
                - supplied.re,
        )?;
        let reactive_power_residual_var =
            finite(supplied.im + (2.0 * self.omega_rad_s) * mean_magnetic_energy_j)?;
        if supplied.re < 0.0
            || power_residual_w.abs() > model.controls.power_tolerance_w
            || reactive_power_residual_var.abs() > model.controls.power_tolerance_w
        {
            return Err(HarmonicWindingError::Balance(
                "terminal power, spatial heat and magnetic storage disagree",
            ));
        }
        poll(cancelled)?;
        Ok(HarmonicWindingResponse {
            omega_rad_s: self.omega_rad_s,
            currents_a: currents_a.to_vec(),
            terminal_voltages_v: voltage,
            potential_wb_m: potential,
            flux_density_t: flux_density,
            eddy_cells: conductors.eddy_cells,
            winding_cells,
            solid_cells: conductors.solid_cells,
            solid_region_currents: conductors.solid_region_currents,
            region_currents_a: conductors.region_currents_a,
            mean_magnetic_energy_j,
            mean_eddy_power_w: conductors.mean_eddy_power_w,
            mean_winding_power_w,
            mean_solid_power_w: conductors.mean_solid_power_w,
            mean_supplied_power_w: supplied.re,
            reactive_power_var: supplied.im,
            power_residual_w,
            reactive_power_residual_var,
        })
    }
}
