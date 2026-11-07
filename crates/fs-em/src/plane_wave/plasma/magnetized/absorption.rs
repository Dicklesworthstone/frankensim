//! Exact finite-volume heat from the original two circular carrier responses.

use super::{PreparedMagnetizedIonizedGasWave, Result, bad, balance, det, positive, sum};

/// Source-bound cycle-mean absorption in a finite prism inside the homogeneous
/// outgoing half-space. This is a control volume, not a film with a second
/// reflecting interface. Signed pool powers already include the original
/// Coulomb drift correction. Zero-drift thermal exchange is a separate process.
pub struct MagnetizedIonizedGasSlabAbsorption<'w, 'a> {
    wave: &'w PreparedMagnetizedIonizedGasWave<'a>,
    area_m2: f64,
    start_depth_m: f64,
    end_depth_m: f64,
    volume_m3: f64,
    electron_power_w: f64,
    heavy_power_w: f64,
    electron_friction_power_w: f64,
    heavy_friction_power_w: f64,
    exchange_correction_to_electrons_w: f64,
    absorbed_power_w: f64,
    drift_kinetic_energy_j: f64,
    background_field_energy_j: f64,
    power_balance_relative_residual: f64,
}

impl<'a> PreparedMagnetizedIonizedGasWave<'a> {
    /// Integrate the actual admitted circular fields through a physical slab.
    /// Requires positive finite area and `0 <= start < end`, with a positive
    /// representable volume. The two analytic attenuation integrals retain
    /// weak absorption and opaque/deep limits without subtracting port powers,
    /// depth quadrature, carrier refactorization or extra field solves.
    ///
    /// The source remains frozen. A finite-energy consumer must bind this
    /// volume to its actual gas parcel and declare the interval over which
    /// these cycle-mean powers are held. Drift kinetic and background field
    /// storage are returned separately and must not be deposited as heat.
    pub fn slab_absorption(
        &self,
        area_m2: f64,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<MagnetizedIonizedGasSlabAbsorption<'_, 'a>> {
        positive(area_m2)?;
        if !start_depth_m.is_finite()
            || !end_depth_m.is_finite()
            || start_depth_m < 0.0
            || end_depth_m <= start_depth_m
        {
            return Err(bad(
                "magnetized absorption needs a finite positive depth interval",
            ));
        }
        let width = end_depth_m - start_depth_m;
        let volume_m3 = positive(area_m2 * width)?;
        let weights = std::array::from_fn::<_, 2, _>(|m| {
            let alpha = -self.interface.modes[m]
                .transmitted_normal_wavenumber_rad_m
                .im;
            SlabWeight::new(volume_m3, alpha, start_depth_m, width)
        });
        let integrate = |value: fn(&super::ModeWave<'a>) -> f64| -> f64 {
            sum((0..2).map(|m| weights[m].integrate(value(&self.modes[m]))))
        };
        let electron_power_w =
            integrate(|m| m.material.mean_heating_w_m3.electron_power_density_w_m3());
        let heavy_power_w = integrate(|m| m.material.mean_heating_w_m3.heavy_power_density_w_m3());
        let electron_friction_power_w = integrate(|m| {
            m.material
                .mean_heating_w_m3
                .electron_friction_power_density_w_m3()
        });
        let heavy_friction_power_w = integrate(|m| {
            m.material
                .mean_heating_w_m3
                .heavy_friction_power_density_w_m3()
        });
        let exchange_correction_to_electrons_w = integrate(|m| {
            m.material
                .mean_heating_w_m3
                .thermal_exchange_correction_to_electrons_w_m3()
        });
        let absorbed_power_w =
            integrate(|m| m.material.mean_heating_w_m3.total_power_density_w_m3());
        let drift_kinetic_energy_j = integrate(|m| m.material.mean_drift_kinetic_energy_j_m3);
        let background_field_energy_j = integrate(|m| m.background_energy_j_m3);
        if [
            electron_power_w,
            heavy_power_w,
            exchange_correction_to_electrons_w,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || [
                electron_friction_power_w,
                heavy_friction_power_w,
                absorbed_power_w,
                drift_kinetic_energy_j,
                background_field_energy_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(bad(
                "magnetized slab heat or stored energy is unrepresentable",
            ));
        }
        let residual = balance(absorbed_power_w, electron_power_w + heavy_power_w)?
            .max(balance(
                absorbed_power_w,
                electron_friction_power_w + heavy_friction_power_w,
            )?)
            .max(self.power_balance_relative_residual);
        Ok(MagnetizedIonizedGasSlabAbsorption {
            wave: self,
            area_m2,
            start_depth_m,
            end_depth_m,
            volume_m3,
            electron_power_w,
            heavy_power_w,
            electron_friction_power_w,
            heavy_friction_power_w,
            exchange_correction_to_electrons_w,
            absorbed_power_w,
            drift_kinetic_energy_j,
            background_field_energy_j,
            power_balance_relative_residual: residual,
        })
    }
}

impl<'w, 'a> MagnetizedIonizedGasSlabAbsorption<'w, 'a> {
    pub const fn wave(&self) -> &'w PreparedMagnetizedIonizedGasWave<'a> {
        self.wave
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn start_depth_m(&self) -> f64 {
        self.start_depth_m
    }
    pub const fn end_depth_m(&self) -> f64 {
        self.end_depth_m
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Signed electron input, including the original Coulomb drift correction.
    pub const fn electron_power_w(&self) -> f64 {
        self.electron_power_w
    }
    /// Signed heavy input; the correction has the opposite sign to electrons.
    pub const fn heavy_power_w(&self) -> f64 {
        self.heavy_power_w
    }
    pub const fn electron_friction_power_w(&self) -> f64 {
        self.electron_friction_power_w
    }
    pub const fn heavy_friction_power_w(&self) -> f64 {
        self.heavy_friction_power_w
    }
    /// Already included in the signed pool powers. Do not deposit twice.
    pub const fn exchange_correction_to_electrons_w(&self) -> f64 {
        self.exchange_correction_to_electrons_w
    }
    pub const fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
    /// Coherent cycle-mean drift storage, separate from deposited heat.
    pub const fn drift_kinetic_energy_j(&self) -> f64 {
        self.drift_kinetic_energy_j
    }
    /// Lossless bound-charge electric and wave magnetic storage, excluding B0.
    pub const fn background_field_energy_j(&self) -> f64 {
        self.background_field_energy_j
    }
    /// Original Maxwell/material residual and integrated thermal-pool closure.
    pub const fn power_balance_relative_residual(&self) -> f64 {
        self.power_balance_relative_residual
    }
}

struct SlabWeight {
    volume: f64,
    average: f64,
    log_average: f64,
    decay: f64,
}
impl SlabWeight {
    fn new(volume: f64, alpha: f64, start: f64, width: f64) -> Self {
        let x = alpha * width;
        let average = if x == 0.0 {
            1.0
        } else {
            (-0.5 * det::expm1(-2.0 * x)) / x
        };
        let log_average = if x.is_infinite() {
            -core::f64::consts::LN_2 - det::ln(alpha) - det::ln(width)
        } else {
            det::ln(average)
        };
        Self {
            volume,
            average,
            log_average,
            decay: -alpha * start,
        }
    }
    fn integrate(&self, value: f64) -> f64 {
        if value == 0.0 {
            return 0.0;
        }
        let a = det::exp(self.decay);
        let direct = (((value * self.volume) * self.average) * a) * a;
        if direct.is_finite() && direct != 0.0 {
            return direct;
        }
        // Rare range fallback: recover a representable physical integral when
        // a large parcel/source and very small attenuation cancel in exponent.
        let log_magnitude =
            det::ln(value.abs()) + det::ln(self.volume) + self.log_average + 2.0 * self.decay;
        det::exp(log_magnitude).copysign(value)
    }
}

#[cfg(test)]
mod tests {
    use super::SlabWeight;

    #[test]
    fn attenuation_integral_retains_zero_and_extreme_exponent_limits() {
        assert_eq!(SlabWeight::new(2.0, 0.0, 1.0, 3.0).integrate(-4.0), -8.0);
        assert_eq!(SlabWeight::new(2.0, 1e-100, 0.0, 3.0).integrate(4.0), 8.0);
        // alpha*width overflows but the complete volume integral is finite.
        let overflow = SlabWeight::new(1e300, 1e200, 0.0, 1e200).integrate(1e200);
        assert!((overflow / 5e99 - 1.0).abs() < 2e-12);
        // Intensity alone underflows while the physical total remains visible.
        let deep = SlabWeight::new(1e200, 1.0, 500.0, 1.0).integrate(1e200);
        let expected =
            (400.0 * core::f64::consts::LN_10 - 1000.0).exp() * (-(-2.0_f64).exp_m1()) / 2.0;
        assert!((deep / expected - 1.0).abs() < 2e-12);
    }
}
