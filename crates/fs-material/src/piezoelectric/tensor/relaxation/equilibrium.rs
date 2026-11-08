//! Static equilibrium of the original declared Maxwell/Debye source model.
//!
//! The auxiliary Debye polarization obeys tau*p_dot+p=delta_epsilon*n.E;
//! its static contribution therefore adds to the instantaneous permittivity.
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_theory.05.36.html>

use super::{RelaxingTensorPiezoelectric, Result, det, invalid, scalar, transient::MemoryResponse};

impl RelaxingTensorPiezoelectric {
    /// Additional equilibrium polarization and memory storage at held field.
    /// Every Maxwell branch has zero relaxed stress and storage. Each Debye
    /// pole contributes p=delta_epsilon*(n.E), polarization n*p, and energy
    /// p²/(2*delta_epsilon). This is STATIC energy, without the harmonic
    /// cycle-average factor of one half; no base-tensor energy is included.
    ///
    /// The original source band must include zero frequency. A finite upper
    /// bound is allowed: this response neither authorizes an unrestricted
    /// transient nor constructs, alters or transfers any transient memory.
    /// The enclosing field owner must solve with the relaxed dielectric
    /// increment and retain these actual equilibrium charges and energies.
    pub fn equilibrium_response(&self, electric_field_v_m: [f64; 3]) -> Result<MemoryResponse> {
        self.admit_frequency(0.0)?;
        if electric_field_v_m.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "piezoelectric equilibrium requires finite electric field",
            ));
        }
        let mut response = MemoryResponse::default();
        for pole in self.dielectric.iter() {
            let along = scalar(
                pole.direction
                    .iter()
                    .zip(electric_field_v_m)
                    .fold(0.0_f64, |sum, (&n, e)| n.mul_add(e, sum)),
            )?;
            // Scale before squaring, as in harmonic memory recovery. Forming
            // E² first would overflow for representable weak-pole energies.
            let root = det::sqrt(pole.permittivity_increment_f_m);
            let normalized = scalar(root * along)?;
            let energy = scalar((0.5 * normalized) * normalized)?;
            let polarization = scalar(root * normalized)?;
            if along != 0.0 && (normalized == 0.0 || energy == 0.0 || polarization == 0.0) {
                return Err(invalid(
                    "piezoelectric equilibrium polarization or storage underflows",
                ));
            }
            response.dielectric_memory_energy_j_m3 =
                scalar(response.dielectric_memory_energy_j_m3 + energy)?;
            for (value, &direction) in response.polarization_c_m2.iter_mut().zip(&pole.direction) {
                *value = scalar(direction.mul_add(polarization, *value))?;
            }
        }
        response.memory_energy_j_m3 = response.dielectric_memory_energy_j_m3;
        Ok(response)
    }
}
