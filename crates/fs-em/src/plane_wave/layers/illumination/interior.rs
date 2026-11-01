//! Local Maxwell fields from the same two-sided finite-film solution.
//!
//! Each parity field uses waves decaying from the nearest physical face. The
//! analytic (1-exp(-i q d))/q limit keeps critical films finite; a bounded
//! reciprocal representation also admits optically deep critical films.
//! No growing exponential, forward/backward division by q, or flux difference
//! is needed for local dielectric heating.

use super::*;

/// Time/ensemble averages at one physical point in a finite film.
/// Axes are x in the incidence plane, y transverse, z from left to right.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LayerPointStatistics {
    /// Componentwise time-averaged E², including the peak-phasor factor 1/2.
    pub electric_rms_squared_v2_m2: [f64; 3],
    /// Componentwise time-averaged H², including the peak-phasor factor 1/2.
    pub magnetic_rms_squared_a2_m2: [f64; 3],
    /// E cross conj(H)/2 [W/m²]: real active and imaginary reactive flux.
    pub complex_poynting_w_m2: [C64; 3],
    /// Positive dielectric loss omega*(-Im epsilon)*<E²> [W/m³].
    /// Includes the longitudinal TM electric field. No thermal state is debited.
    pub absorbed_power_w_m3: f64,
}

/// Peak coherent Maxwell fields and their local time-averaged observables.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerPointField {
    pub electric_v_m: [C64; 3],
    pub magnetic_a_m: [C64; 3],
    pub statistics: LayerPointStatistics,
}

/// Two source columns at one authored film/depth, bound to the original
/// illumination and its material, frequency, transverse wavevector and clock.
/// Preparation and coherent evaluation cost O(1); neither allocates. Repeated
/// coherence evaluation reuses these fields without repeating any trigonometry.
pub struct PreparedLayerPoint<'p, 's, 'm> {
    illumination: &'p PreparedLayerIllumination<'s, 'm>,
    layer_index: usize,
    depth_m: f64,
    distance_from_right_m: f64,
    electric: [[C64; 2]; 3],
    magnetic: [[C64; 2]; 3],
    root_loss: f64,
}

// Return p and 1-p without cancelling weak attenuation or a small phase.
// Fully decayed waves do not need an otherwise unrepresentable optical phase.
fn propagation(q: C64, distance: f64) -> Result<(C64, C64), ReflectionError> {
    if !distance.is_finite() || distance < 0.0 {
        return Err(bad("interior optical distance is unrepresentable"));
    }
    let decay = q.im * distance;
    let attenuation = det::exp(decay);
    if attenuation == 0.0 {
        return Ok((C64::ZERO, C64::ONE));
    }
    let phase = q.re * distance;
    if !phase.is_finite() || !attenuation.is_finite() || attenuation > 1.0 {
        return Err(bad("interior propagation is unrepresentable"));
    }
    let (sine, cosine) = det::sin_cos(phase);
    let half_sine = det::sin(0.5 * phase);
    Ok((
        C64::new(attenuation * cosine, -attenuation * sine),
        C64::new(
            -det::expm1(decay) + 2.0 * attenuation * half_sine * half_sine,
            attenuation * sine,
        ),
    ))
}

impl<'s, 'm> PreparedLayerIllumination<'s, 'm> {
    /// Prepare an interior point, measured from this film's left physical face.
    /// Face traces at 0 and thickness are allowed on the selected film side.
    /// Zero-thickness identities have no interior volume and refuse this query.
    pub fn prepare_point(
        &self,
        layer_index: usize,
        depth_m: f64,
    ) -> Result<PreparedLayerPoint<'_, 's, 'm>, ReflectionError> {
        self.point_from_face(layer_index, depth_m, false)
    }

    /// Reference a point to the right film face. This preserves shallow exit
    /// depths even when subtracting them from an enormous thickness would round
    /// back to that thickness. The physical x/y/z axes are unchanged.
    pub fn prepare_point_from_right(
        &self,
        layer_index: usize,
        distance_m: f64,
    ) -> Result<PreparedLayerPoint<'_, 's, 'm>, ReflectionError> {
        self.point_from_face(layer_index, distance_m, true)
    }

    fn point_from_face(
        &self,
        layer_index: usize,
        distance_m: f64,
        from_right: bool,
    ) -> Result<PreparedLayerPoint<'_, 's, 'm>, ReflectionError> {
        let layer = self
            .stack
            .layers
            .get(layer_index)
            .ok_or_else(|| bad("interior film index is out of range"))?;
        let thickness = layer.source.thickness_m;
        if thickness == 0.0 || !distance_m.is_finite() || distance_m < 0.0 || distance_m > thickness
        {
            return Err(bad(
                "interior depth needs a positive film and lies outside it",
            ));
        }
        let (depth_m, distance_from_right_m) = if from_right {
            (thickness - distance_m, distance_m)
        } else {
            (distance_m, thickness - distance_m)
        };
        let c = self.incidence_cosine;
        let q = layer.wave.normal(c)?;
        let normal_squared = layer.wave.normal_squared(c);
        // dE/d(k_inc z)=-i a H/Y_left, d(H/Y_left)/d(k_inc z)=-i b E.
        let (a, b) = match self.polarization {
            Polarization::Te => {
                let a = C64::from_re(layer.wave.mu_ratio * c);
                (a, normal_squared / a)
            }
            Polarization::Tm => {
                let b = layer.wave.epsilon_ratio.scale(c);
                (normal_squared / b, b)
            }
        };
        if !finite(a) || !finite(b) {
            return Err(bad("interior Maxwell coupling is unrepresentable"));
        }
        let k = self.stack.wave_number;
        let distance = positive(k * thickness)?;
        let (p, difference) = propagation(q, distance)?;
        let (left, _) = propagation(q, k * depth_m)?;
        let (right, _) = propagation(q, k * distance_from_right_m)?;
        // [u,v] represents [1,(1-p)/q] up to a common scale, without
        // overflowing at a deep critical film or dividing a zero normal wave.
        let reciprocal = if q == C64::ZERO {
            distance > 1.0
        } else {
            maximum(difference) > maximum(q)
        };
        let (u, v) = if q == C64::ZERO {
            if reciprocal {
                (C64::new(0.0, -1.0 / distance), C64::ONE)
            } else {
                (C64::ONE, C64::new(0.0, distance))
            }
        } else if reciprocal {
            (q / difference, C64::ONE)
        } else {
            (C64::ONE, difference / q)
        };
        // u*(left-right)/q. Use the nearest face, so an opaque film never
        // subtracts two enormous midpoint phases to recover a shallow depth.
        let local_difference = if q == C64::ZERO {
            if reciprocal {
                C64::from_re((distance_from_right_m - depth_m) / thickness)
            } else {
                C64::new(0.0, k * (distance_from_right_m - depth_m))
            }
        } else {
            let (near, span, sign) = if depth_m <= distance_from_right_m {
                (left, distance_from_right_m - depth_m, 1.0)
            } else {
                (right, depth_m - distance_from_right_m, -1.0)
            };
            let (_, local) = propagation(q, k * span)?;
            near * (local / if reciprocal { difference } else { q }).scale(sign)
        };
        let sum = left + right;
        let parity = |coupling: C64| {
            let scale = maximum(coupling).max(1.0);
            let reference = divided(u, scale);
            let coupling = divided(coupling, scale);
            let denominator = (C64::ONE + p) * reference + coupling * v;
            [
                sum * reference / denominator,
                coupling * local_difference / denominator,
            ]
        };
        let [e_even, h_even] = parity(b);
        let [h_odd, e_odd] = parity(a);
        let inputs = self.interior_inputs[layer_index];
        let mut electric = [[C64::ZERO; 2]; 3];
        let mut magnetic = [[C64::ZERO; 2]; 3];
        let sine = det::sqrt((1.0 - c) * (1.0 + c));
        let base_admittance = match self.polarization {
            Polarization::Te => self.admittance_s[0] / c,
            Polarization::Tm => self.admittance_s[0] * c,
        };
        for source in 0..2 {
            let even = inputs[0][source] + inputs[1][source];
            let odd = inputs[0][source] - inputs[1][source];
            let e = dot([e_even, e_odd], [even, odd]);
            let h = dot([h_even, h_odd], [even, odd]).scale(self.admittance_s[0]);
            match self.polarization {
                Polarization::Te => {
                    electric[1][source] = e;
                    magnetic[0][source] = -h;
                    magnetic[2][source] = e.scale(base_admittance * (sine / layer.wave.mu_ratio));
                }
                Polarization::Tm => {
                    electric[0][source] = e;
                    electric[2][source] =
                        (C64::from_re(-sine / base_admittance) / layer.wave.epsilon_ratio) * h;
                    magnetic[1][source] = h;
                }
            }
        }
        let root_loss = det::sqrt(self.stack.omega)
            * det::sqrt(-layer.wave.epsilon_f_m.im)
            * core::f64::consts::FRAC_1_SQRT_2;
        if !root_loss.is_finite()
            || electric
                .iter()
                .chain(&magnetic)
                .flatten()
                .any(|&z| !finite(z))
        {
            return Err(bad("interior field transfer is unrepresentable"));
        }
        Ok(PreparedLayerPoint {
            illumination: self,
            layer_index,
            depth_m,
            distance_from_right_m,
            electric,
            magnetic,
            root_loss,
        })
    }
}

impl<'p, 's, 'm> PreparedLayerPoint<'p, 's, 'm> {
    pub fn illumination(&self) -> &'p PreparedLayerIllumination<'s, 'm> {
        self.illumination
    }
    pub fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub fn depth_m(&self) -> f64 {
        self.depth_m
    }
    pub fn distance_from_right_m(&self) -> f64 {
        self.distance_from_right_m
    }

    /// Contract the original simultaneous peak tangential E sources [V/m].
    /// Standing-wave cancellation happens before any square or heat evaluation.
    pub fn evaluate(&self, input: [C64; 2]) -> Result<LayerPointField, ReflectionError> {
        if input.iter().any(|&z| !finite(z)) {
            return Err(bad("nonfinite interior illumination source"));
        }
        let electric = self.electric.map(|row| dot(row, input));
        let magnetic = self.magnetic.map(|row| dot(row, input));
        let rms_squared =
            |field: [C64; 3]| field.map(|z| z.scale(core::f64::consts::FRAC_1_SQRT_2).norm_sq());
        let poynting = std::array::from_fn(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            dot(
                [electric[j], -electric[k]],
                [magnetic[k].conj(), magnetic[j].conj()],
            )
            .scale(0.5)
        });
        let heat = electric
            .iter()
            .fold(Dd::ZERO, |s, &e| {
                s + Dd::from_f64(e.scale(self.root_loss).norm_sq())
            })
            .to_f64();
        let statistics = LayerPointStatistics {
            electric_rms_squared_v2_m2: rms_squared(electric),
            magnetic_rms_squared_a2_m2: rms_squared(magnetic),
            complex_poynting_w_m2: poynting,
            absorbed_power_w_m3: heat,
        };
        if electric
            .iter()
            .chain(&magnetic)
            .chain(&poynting)
            .any(|&z| !finite(z))
            || statistics
                .electric_rms_squared_v2_m2
                .iter()
                .chain(&statistics.magnetic_rms_squared_a2_m2)
                .chain(std::iter::once(&heat))
                .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(bad("interior field or power overflow"));
        }
        Ok(LayerPointField {
            electric_v_m: electric,
            magnetic_a_m: magnetic,
            statistics,
        })
    }

    /// Explicit peak-E second-moment factors C=B B*, each in V/m. Independent
    /// factors add observables, never phasors. These are spectral-line averages,
    /// not PSD/Hz. This point counts as one layer in maximum_layer_mode_pairs.
    /// Allocation-free and transactional; zero factors mean exact silence.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[[C64; 2]],
        budget: IlluminationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayerPointStatistics, ReflectionError> {
        if factors.len() > 4096
            || factors.len() > budget.maximum_modes
            || factors.len() > budget.maximum_layer_mode_pairs
        {
            return Err(bad("interior coherence exceeds work budget"));
        }
        if cancelled() {
            return Err(bad("illumination cancelled"));
        }
        let mut sums = [Dd::ZERO; 13];
        for &factor in factors {
            if cancelled() {
                return Err(bad("illumination cancelled"));
            }
            let value = self.evaluate(factor)?.statistics;
            for i in 0..3 {
                sums[i] = sums[i] + Dd::from_f64(value.electric_rms_squared_v2_m2[i]);
                sums[3 + i] = sums[3 + i] + Dd::from_f64(value.magnetic_rms_squared_a2_m2[i]);
                sums[6 + 2 * i] = sums[6 + 2 * i] + Dd::from_f64(value.complex_poynting_w_m2[i].re);
                sums[7 + 2 * i] = sums[7 + 2 * i] + Dd::from_f64(value.complex_poynting_w_m2[i].im);
            }
            sums[12] = sums[12] + Dd::from_f64(value.absorbed_power_w_m3);
        }
        let total = sums.map(Dd::to_f64);
        if total.iter().any(|x| !x.is_finite()) {
            return Err(bad("interior statistical power overflow"));
        }
        Ok(LayerPointStatistics {
            electric_rms_squared_v2_m2: std::array::from_fn(|i| total[i]),
            magnetic_rms_squared_a2_m2: std::array::from_fn(|i| total[3 + i]),
            complex_poynting_w_m2: std::array::from_fn(|i| {
                C64::new(total[6 + 2 * i], total[7 + 2 * i])
            }),
            absorbed_power_w_m3: total[12],
        })
    }

    /// Independent face illumination is an explicit diagonal source covariance.
    pub fn evaluate_independent_sources(
        &self,
        input: [C64; 2],
    ) -> Result<LayerPointStatistics, ReflectionError> {
        self.evaluate_coherence_modes(
            &[[input[0], C64::ZERO], [C64::ZERO, input[1]]],
            IlluminationBudget {
                maximum_modes: 2,
                maximum_layer_mode_pairs: 2,
            },
            || false,
        )
    }
}
