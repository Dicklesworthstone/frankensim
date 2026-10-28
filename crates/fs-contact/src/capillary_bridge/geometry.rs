//! Exact circular Young--Laplace profiles in the declared two-dimensional lane.
use super::{BridgeError, BridgeProfile, finite, positive};
use fs_material::liquid::surface::adhesion::YoungDupreLaw;
use fs_math::det;

#[derive(Clone, Copy, Debug)]
pub(super) struct Geometry {
    pub volume: f64,
    pub width: f64,
    pub area: f64,
    pub sigma_width: f64,
    pub theta: f64,
    pub cosine: f64,
    pub sine: f64,
    pub alpha: f64,
    pub segment: f64,
    pub neck: f64,
    pub energy_slope: f64,
    pub rupture_gap: f64,
    pub formation_gap: f64,
    pub cap_radius: f64,
    pub cap_energy: f64,
}

impl Geometry {
    pub fn new(volume: f64, width: f64, law: YoungDupreLaw) -> Result<Self, BridgeError> {
        positive(volume, "liquid volume")?;
        positive(width, "invariant extrusion width")?;
        if volume * 0.5 <= 0.0 || (volume * 0.5) * 2.0 != volume {
            return Err(BridgeError::Invalid {
                field: "two representable equal cap volumes",
            });
        }
        let area = positive(volume / width, "liquid cross-sectional area")?;
        let sigma = law.surface_tension_n_m();
        // Use the actual thermodynamic coefficient so geometry and force share
        // the SAME wetted-surface energy, also for an adhesion-derived angle.
        let cosine = -law.wetted_surface_energy_difference_j_m2() / sigma;
        if !(cosine > 0.0 && cosine < 1.0) {
            return Err(BridgeError::Invalid {
                field: "equal acute partial-wetting angles",
            });
        }
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
        let theta = det::atan2(sine, cosine);
        let alpha = det::atan2(cosine, sine);
        let segment = angle_minus_sin_cos(alpha) / (2.0 * cosine * cosine);
        let neck = neck_numerator(alpha) / (2.0 * cosine * cosine);
        let energy_slope = alpha / cosine + sine;
        positive(segment, "meniscus segment coefficient")?;
        positive(neck, "meniscus neck coefficient")?;
        let rupture_gap = positive(det::sqrt(area / neck), "geometric rupture gap")?;
        let cap_segment = positive(angle_minus_sin_cos(theta), "detached cap segment")?;
        let cap_radius = positive(det::sqrt((0.5 * area) / cap_segment), "detached cap radius")?;
        let half_sine = det::sin(0.5 * theta);
        let formation_gap = positive(4.0 * cap_radius * half_sine * half_sine, "cap touching gap")?;
        let sigma_width = positive(sigma * width, "surface tension times extrusion width")?;
        let cap_energy = positive(
            4.0 * sigma_width * cap_radius * cap_segment,
            "two-cap surface energy",
        )?;
        if formation_gap >= rupture_gap {
            return Err(BridgeError::Invalid {
                field: "ordered formation and geometric rupture gaps",
            });
        }
        Ok(Self {
            volume,
            width,
            area,
            sigma_width,
            theta,
            cosine,
            sine,
            alpha,
            segment,
            neck,
            energy_slope,
            rupture_gap,
            formation_gap,
            cap_radius,
            cap_energy,
        })
    }

    pub fn attached(self, gap: f64) -> Result<BridgeProfile, BridgeError> {
        positive(gap, "plate gap")?;
        if gap > self.rupture_gap {
            return Err(BridgeError::NoConnectedProfile);
        }
        let span = finite(self.area / gap + self.segment * gap, "wetted span")?;
        // Factored at the known geometric endpoint; avoids cancellation of two
        // almost equal lengths and represents exactly zero at the topology rule.
        let neck = finite(
            self.neck * (self.rupture_gap - gap) * ((self.rupture_gap + gap) / gap),
            "bridge neck",
        )?;
        let force = finite(
            self.sigma_width * (2.0 * self.cosine * (self.area / gap) / gap + self.energy_slope),
            "bridge force",
        )?;
        Ok(BridgeProfile {
            gap_m: gap,
            attached: true,
            liquid_volume_m3: self.volume,
            detached_cap_volumes_m3: None,
            wetted_span_m: span,
            neck_width_m: Some(neck),
            liquid_gas_area_m2: finite(
                2.0 * gap * (self.alpha / self.cosine) * self.width,
                "liquid-gas area",
            )?,
            wetted_area_m2: finite(2.0 * span * self.width, "wetted area")?,
            liquid_minus_gas_pressure_pa: finite(
                -2.0 * (self.sigma_width / self.width) * self.cosine / gap,
                "Laplace pressure",
            )?,
            surface_free_energy_j: finite(
                self.sigma_width
                    * (-2.0 * self.cosine * (self.area / gap) + self.energy_slope * gap),
                "bridge surface energy",
            )?,
            normal_force_n: -force,
            force_gap_derivative_n_per_m: finite(
                4.0 * self.sigma_width * self.cosine * (self.area / gap) / gap / gap,
                "bridge force derivative",
            )?,
        })
    }

    pub fn detached(self, gap: f64) -> Result<BridgeProfile, BridgeError> {
        positive(gap, "plate gap")?;
        if gap < self.formation_gap {
            return Err(BridgeError::OverlappingCaps);
        }
        let span = finite(2.0 * self.cap_radius * self.sine, "detached wetted span")?;
        Ok(BridgeProfile {
            gap_m: gap,
            attached: false,
            liquid_volume_m3: self.volume,
            detached_cap_volumes_m3: Some([0.5 * self.volume; 2]),
            wetted_span_m: span,
            neck_width_m: None,
            liquid_gas_area_m2: finite(
                4.0 * self.cap_radius * self.theta * self.width,
                "detached liquid-gas area",
            )?,
            wetted_area_m2: finite(2.0 * span * self.width, "detached wetted area")?,
            liquid_minus_gas_pressure_pa: finite(
                (self.sigma_width / self.width) / self.cap_radius,
                "detached Laplace pressure",
            )?,
            surface_free_energy_j: self.cap_energy,
            normal_force_n: 0.0,
            force_gap_derivative_n_per_m: 0.0,
        })
    }

    /// Exact fixed-volume work along the connected branch, avoiding subtraction
    /// of nearly equal free energies. Positive work increases plate separation.
    pub fn work(self, from: f64, to: f64) -> Result<f64, BridgeError> {
        let work = finite(
            (to - from)
                * (self.sigma_width
                    * (2.0 * self.cosine * (self.area / from) / to + self.energy_slope)),
            "capillary work",
        )?;
        if from != to && work == 0.0 {
            return Err(BridgeError::Invalid {
                field: "unrepresentable nonzero capillary work",
            });
        }
        Ok(work)
    }
}

fn angle_minus_sin_cos(x: f64) -> f64 {
    if x < 0.01 {
        let x2 = x * x;
        x * x2
            * (2.0 / 3.0
                + x2 * (-2.0 / 15.0
                    + x2 * (4.0 / 315.0 + x2 * (-2.0 / 2835.0 + x2 * 4.0 / 155925.0))))
    } else {
        x - det::sin(x) * det::cos(x)
    }
}
fn neck_numerator(x: f64) -> f64 {
    if x < 0.01 {
        let x2 = x * x;
        x * x2
            * (1.0 / 3.0
                + x2 * (-7.0 / 60.0
                    + x2 * (31.0 / 2520.0 + x2 * (-127.0 / 181440.0 + x2 * 511.0 / 19958400.0))))
    } else {
        2.0 * det::sin(x) - det::sin(x) * det::cos(x) - x
    }
}
