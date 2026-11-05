//! Gray surface exchange through the original two source-integrated layers.

use super::*;

// Same rounded CODATA value as the original thermal gray-surface owners.
const SIGMA: f64 = 5.670_374_419e-8;

/// Two opaque diffuse-gray, equal-area faces with full mutual view, separated
/// by a nonparticipating gap. Gas/contact conduction acts in parallel with
/// radiation BETWEEN these solved surfaces, not between the bulk cell states.
/// The fixed HTC excludes both material half cells and all radiative transfer.
/// No finite gap storage, view-factor inference or emissivity evolution occurs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayLayerGap {
    effective_emissivity: f64,
    conductance_w_m2_k: f64,
}

impl GrayLayerGap {
    /// Explicitly declare F_A_to_B = F_B_to_A = 1 for equal-area faces.
    /// Emissivities alone do not establish this geometric/view assumption.
    pub fn equal_area_full_view(
        emissivities: [f64; 2],
        conductance_w_m2_k: f64,
    ) -> Result<Self, DepositionError> {
        if emissivities
            .iter()
            .any(|e| !e.is_finite() || !(0.0..=1.0).contains(e))
            || !conductance_w_m2_k.is_finite()
            || conductance_w_m2_k < 0.0
        {
            return Err(invalid(
                "gray layer gap needs emissivities in [0,1] and nonnegative finite HTC",
            ));
        }
        let low = emissivities[0].min(emissivities[1]);
        let high = emissivities[0].max(emissivities[1]);
        let effective_emissivity = if low == 0.0 {
            0.0
        } else {
            low / ((1.0 - low) + low / high)
        };
        Ok(Self {
            effective_emissivity,
            conductance_w_m2_k,
        })
    }

    fn at(self, surface_k: [f64; 2]) -> Result<GapState, DepositionError> {
        let [a, b] = surface_k;
        let delta = a - b;
        let radiative = if self.effective_emissivity == 0.0 {
            0.0
        } else {
            positive(SIGMA * self.effective_emissivity)?
        };
        let radiation = if radiative == 0.0 || delta == 0.0 {
            0.0
        } else {
            // Factor the fourth-power difference to retain near equilibrium
            // heat without subtracting two large absolute radiant powers.
            finite(((radiative * delta) * (a + b)) * finite(a * a + b * b)?)?
        };
        let conduction = finite(self.conductance_w_m2_k * delta)?;
        let [a_slope, b_slope] =
            surface_k.map(|t| finite(self.conductance_w_m2_k + ((4.0 * radiative * t) * t) * t));
        Ok(GapState {
            conduction_w_m2: conduction,
            radiation_w_m2: radiation,
            slopes_w_m2_k: [a_slope?, b_slope?],
        })
    }
}

struct GapState {
    conduction_w_m2: f64,
    radiation_w_m2: f64,
    // Positive magnitudes: d(j)/dTa and -d(j)/dTb.
    slopes_w_m2_k: [f64; 2],
}

/// One conservative flux, positive from original A to B. Surface powers are
/// two components of that SAME transfer and must not be credited separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayLayerConductionResponse {
    pub power_a_to_b_w: f64,
    pub conductance_w_k: f64,
    /// Ordered material-face temperatures after the source-integrated drops.
    pub interface_temperature_k: [f64; 2],
    /// Full derivatives [dQ/dTa, dQ/dTb], with material charts and gap fixed.
    pub endpoint_derivative_w_k: [f64; 2],
    pub gap_conduction_power_w: f64,
    pub radiation_power_w: f64,
    /// Hot-layer outgoing minus cold-layer incoming power [W].
    pub layer_heat_residual_w: f64,
    /// Original A-to-B power minus both gap contributions [W]. When radiation
    /// is disabled, the original contact's Kelvin/roundoff gate remains
    /// authoritative near perfect contact rather than this ill-conditioned
    /// diagnostic formed from a very small surface temperature jump.
    pub surface_heat_residual_w: f64,
    pub iterations: usize,
}

impl TwoLayerConduction<'_> {
    /// Integrate each original k(T) profile to its solved gray-gap surface.
    ///
    /// The existing scalar flux bracket, exact piecewise-affine Kirchhoff
    /// inverses and safeguarded Newton iteration are shared with `evaluate`.
    /// The interface equation is j = h_g*(Ua-Ub) + epsilon_eff*sigma*(Ua^4-Ub^4).
    /// Source support is required only inside each own material half cell;
    /// disjoint temperature charts may be bridged by the physical gap.
    ///
    /// Authored contact resistance must be zero: the gap owns the entire
    /// interface. Bulk layers remain source-free, quasi-steady pure phases;
    /// no melting profile, surface storage, coating or extra heat is inferred.
    /// A consuming thermal solver still owns both finite bodies and the clock.
    pub fn evaluate_gray_gap(
        &self,
        temperatures_k: [f64; 2],
        gap: GrayLayerGap,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayLayerConductionResponse, DepositionError> {
        checkpoint()?;
        self.admit_evaluation(temperatures_k, config)?;
        if self.contact_resistance_m2_k_w != 0.0 {
            return Err(invalid(
                "gray layer gap must replace zero authored contact resistance",
            ));
        }
        let [ta, tb] = temperatures_k;
        let disconnected = gap.effective_emissivity == 0.0 && gap.conductance_w_m2_k == 0.0;
        if gap.effective_emissivity == 0.0 && !disconnected {
            // Preserve the exact original contact result and its near-perfect
            // contact representation gate, not merely its mathematical limit.
            let model = Self {
                contact_resistance_m2_k_w: positive(gap.conductance_w_m2_k.recip())?,
                ..*self
            };
            let response = model.evaluate(temperatures_k, config, &mut checkpoint)?;
            let gap_power = finite(response.power_a_to_b_w + response.contact_heat_residual_w)?;
            return Ok(GrayLayerConductionResponse {
                power_a_to_b_w: response.power_a_to_b_w,
                conductance_w_k: response.conductance_w_k,
                interface_temperature_k: response.interface_temperature_k,
                endpoint_derivative_w_k: response.endpoint_derivative_w_k,
                gap_conduction_power_w: gap_power,
                radiation_power_w: 0.0,
                layer_heat_residual_w: response.layer_heat_residual_w,
                surface_heat_residual_w: finite(response.power_a_to_b_w - gap_power)?,
                iterations: response.iterations,
            });
        }
        if disconnected || ta == tb {
            let g = if disconnected {
                0.0
            } else {
                let transfer = positive(gap.at(temperatures_k)?.slopes_w_m2_k[0])?;
                series_contact_conductance_w_k(
                    self.area_m2,
                    self.distances_m,
                    [
                        self.curves[0].conductivity_w_m_k(ta)?,
                        self.curves[1].conductivity_w_m_k(tb)?,
                    ],
                    positive(transfer.recip())?,
                )?
            };
            checkpoint()?;
            return Ok(GrayLayerConductionResponse {
                power_a_to_b_w: 0.0,
                conductance_w_k: g,
                interface_temperature_k: temperatures_k,
                endpoint_derivative_w_k: [g, -g],
                gap_conduction_power_w: 0.0,
                radiation_power_w: 0.0,
                layer_heat_residual_w: 0.0,
                surface_heat_residual_w: 0.0,
                iterations: 0,
            });
        }
        let (hot, cold, sign) = if ta > tb { (0, 1, 1.0) } else { (1, 0, -1.0) };
        let (th, tc) = (temperatures_k[hot], temperatures_k[cold]);
        let (kh, kc) = (self.curves[hot], self.curves[cold]);
        let (dh, dc) = (self.distances_m[hot], self.distances_m[cold]);
        self.solve_interface(temperatures_k, config, None, &mut checkpoint, |trial| {
            let [u, v] = trial.surface_temperature_k;
            let [jh, jc] = trial.layer_flux_w_m2;
            let surface = gap.at([u, v])?;
            let interface_flux = finite(surface.conduction_w_m2 + surface.radiation_w_m2)?;
            let qh = finite(self.area_m2 * jh)?;
            let qc = finite(self.area_m2 * jc)?;
            let gap_power = finite(self.area_m2 * surface.conduction_w_m2)?;
            let radiation_power = finite(self.area_m2 * surface.radiation_w_m2)?;
            let surface_power = finite(gap_power + radiation_power)?;
            let layer_residual = finite(qh - qc)?;
            let surface_residual = finite(qh - surface_power)?;
            let tolerance = finite(
                config.absolute_heat_tolerance_w
                    + config.relative_heat_tolerance
                        * qh.abs().max(qc.abs()).max(surface_power.abs()),
            )?;
            let ku = kh.conductivity_w_m_k(u)?;
            let kv = kc.conductivity_w_m_k(v)?;
            let [a, b] = surface.slopes_w_m2_k;
            let denominator = positive(1.0 + a * (dh / ku) + b * (dc / kv))?;
            if layer_residual.abs() <= tolerance && surface_residual.abs() <= tolerance {
                let mut interface = [0.0; 2];
                interface[hot] = u;
                interface[cold] = v;
                let mut derivatives = [0.0; 2];
                derivatives[hot] = sign
                    * positive(
                        (self.area_m2 * (a / denominator)) * (kh.conductivity_w_m_k(th)? / ku),
                    )?;
                derivatives[cold] = -sign
                    * positive(
                        (self.area_m2 * (b / denominator)) * (kc.conductivity_w_m_k(tc)? / kv),
                    )?;
                return Ok(LayerFluxUpdate::Accepted(GrayLayerConductionResponse {
                    power_a_to_b_w: sign * qh,
                    conductance_w_k: positive(qh / (th - tc))?,
                    interface_temperature_k: interface,
                    endpoint_derivative_w_k: derivatives,
                    gap_conduction_power_w: sign * gap_power,
                    radiation_power_w: sign * radiation_power,
                    layer_heat_residual_w: layer_residual,
                    surface_heat_residual_w: sign * surface_residual,
                    iterations: trial.iterations,
                }));
            }
            Ok(LayerFluxUpdate::Continue {
                decreasing_residual: finite(interface_flux - trial.flux_w_m2)?,
                slope_magnitude: denominator,
            })
        })
    }
}
