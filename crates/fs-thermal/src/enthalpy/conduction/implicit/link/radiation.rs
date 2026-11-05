//! Two finite bodies exchanging heat through one pair of gray gap surfaces.

use super::super::{ImplicitConductionLaw, ImplicitConductionState};
use super::{ImplicitLinkLaw, ImplicitLinkState};
use crate::enthalpy::conduction::ConductiveLink;
use crate::enthalpy::{ContactPatch, DepositionError, ThermalLaw};

// Same rounded CODATA value as the original fs-conduction radiation owner.
const SIGMA: f64 = 5.670_374_419e-8;

/// Bounded surface-balance work, independent of thermal timestep accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayRadiativeControls {
    /// Maximum safeguarded scalar evaluations, in 1..=128.
    pub maximum_iterations: usize,
    /// Absolute original surface power-balance allowance [W].
    pub absolute_power_tolerance_w: f64,
    /// Relative allowance on the actual exchanged powers, in [0,1).
    pub relative_power_tolerance: f64,
}
impl Default for GrayRadiativeControls {
    fn default() -> Self {
        Self {
            maximum_iterations: 64,
            absolute_power_tolerance_w: 1e-10,
            relative_power_tolerance: 1e-12,
        }
    }
}

/// One complete face solution, with positive power from original A to B.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayRadiativeResponse {
    /// Native conservative secant and complete ordered enthalpy derivatives.
    pub link: ImplicitLinkState,
    /// Ordered exterior temperatures [K], after the two half-cell/coating drops.
    /// None means an insulating bulk leaves a massless face undetermined.
    pub surface_temperatures_k: Option<[f64; 2]>,
    /// Exactly the power applied by the native link, G * (T_A - T_B) [W].
    pub power_w: f64,
    /// Conductive gas-gap contribution between the SAME two surfaces [W].
    pub gap_power_w: f64,
    /// Net gray radiation between those surfaces [W].
    pub radiation_power_w: f64,
    /// Applied power minus the sum of both surface contributions [W].
    pub surface_residual_w: f64,
    /// Actual scalar evaluations; zero for an insulating/disconnected link.
    pub iterations: usize,
}

/// Conservative radiation and parallel gas-gap conduction between finite cells.
///
/// The two opaque, diffuse-gray faces have equal area and full mutual view;
/// radiation does not attenuate in the gap. The two-surface radiosity equations
/// reduce to epsilon_eff = 1 / (1/epsilon_A + 1/epsilon_B - 1). See the original
/// theory at <https://mooseframework.inl.gov/modules/heat_transfer/>.
/// Fixed emissivities, gap HTC and two coating resistances are explicit inputs.
///
/// Each original half-cell retains its own cellwise k(h), including its complete
/// enthalpy derivative on latent plateaus. Its coating is in series BEFORE its
/// radiating surface. The geometry's unsplit contact resistance must be zero;
/// it cannot be assigned to either side silently. The nominal and both covered
/// patch areas must be equal physical areas. This law replaces the selected
/// link through `with_link_laws`; it must not also be installed as a source.
///
/// Both bodies keep their original finite thermal storage. The existing link
/// owner applies one signed power with equal/opposite cell contributions and
/// its original entropy production. No prescribed bath, extra surface storage,
/// enclosure solver, profile-integrated k(T), source-dependent emissivity,
/// near-field radiation or participating-medium transport is introduced.
/// Every successful evaluation is allocation-free.
pub struct GrayRadiativeLink<'a> {
    materials: [&'a dyn ImplicitConductionLaw; 2],
    effective_emissivity: f64,
    gap_htc_w_m2_k: f64,
    coating_resistances_m2_k_w: [f64; 2],
    controls: GrayRadiativeControls,
}

impl<'a> GrayRadiativeLink<'a> {
    /// Bind the SAME ordered caloric/conductivity objects as the finite cells.
    pub fn new(
        materials: [&'a dyn ImplicitConductionLaw; 2],
        emissivities: [f64; 2],
        gap_htc_w_m2_k: f64,
        coating_resistances_m2_k_w: [f64; 2],
        controls: GrayRadiativeControls,
    ) -> Result<Self, DepositionError> {
        if emissivities
            .iter()
            .any(|e| !e.is_finite() || !(0.0..=1.0).contains(e))
            || !nonnegative(gap_htc_w_m2_k)
            || coating_resistances_m2_k_w.iter().any(|r| !nonnegative(*r))
            || !(1..=128).contains(&controls.maximum_iterations)
            || !nonnegative(controls.absolute_power_tolerance_w)
            || !nonnegative(controls.relative_power_tolerance)
            || controls.relative_power_tolerance >= 1.0
        {
            return Err(DepositionError::InvalidInput(
                "gray link needs emissivities in [0,1], nonnegative gap/coatings and bounded controls",
            ));
        }
        // Do not form epsilon_A * epsilon_B or reciprocal tiny emissivities.
        let low = emissivities[0].min(emissivities[1]);
        let high = emissivities[0].max(emissivities[1]);
        let effective_emissivity = if low == 0.0 {
            0.0
        } else {
            low / ((1.0 - low) + low / high)
        };
        Ok(Self {
            materials,
            effective_emissivity,
            gap_htc_w_m2_k,
            coating_resistances_m2_k_w,
            controls,
        })
    }

    /// Resolve both surfaces, their original balance and analytic tangents.
    /// Refusal or cancellation returns no partial face solution.
    pub fn evaluate(
        &self,
        enthalpies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayRadiativeResponse, DepositionError> {
        checkpoint()?;
        let area = patch.link.area_m2;
        if geometry.link_id != patch.link.id
            || patch.link.a == patch.link.b
            || !positive(area)
            || patch.covered_area_m2 != [area; 2]
            || geometry.half_distances_m.iter().any(|d| !positive(*d))
            || geometry.contact_resistance_m2_k_w != 0.0
        {
            return Err(DepositionError::InvalidInput(
                "gray link needs equal physical face areas and explicit per-side coating resistance",
            ));
        }
        let states = [
            self.state(0, enthalpies_j_kg[0], checkpoint)?,
            self.state(1, enthalpies_j_kg[1], checkpoint)?,
        ];
        let temperatures = [states[0].temperature_k, states[1].temperature_k];
        let conductivities = [states[0].conductivity_w_m_k, states[1].conductivity_w_m_k];
        let disconnected = self.effective_emissivity == 0.0 && self.gap_htc_w_m2_k == 0.0;
        if disconnected || conductivities.contains(&0.0) {
            let surfaces = match (conductivities[0] == 0.0, conductivities[1] == 0.0) {
                (false, false) => Some(temperatures),
                (true, false) if !disconnected => Some([temperatures[1]; 2]),
                (false, true) if !disconnected => Some([temperatures[0]; 2]),
                _ => None,
            };
            checkpoint()?;
            return Ok(GrayRadiativeResponse {
                link: ImplicitLinkState {
                    conductance_w_k: 0.0,
                    power_enthalpy_derivatives_w_kg_j: [0.0; 2],
                },
                surface_temperatures_k: surfaces,
                power_w: 0.0,
                gap_power_w: 0.0,
                radiation_power_w: 0.0,
                surface_residual_w: 0.0,
                iterations: 0,
            });
        }
        let bulk_resistances = [
            finite(geometry.half_distances_m[0] / conductivities[0])?,
            finite(geometry.half_distances_m[1] / conductivities[1])?,
        ];
        let resistances = [
            finite(bulk_resistances[0] + self.coating_resistances_m2_k_w[0])?,
            finite(bulk_resistances[1] + self.coating_resistances_m2_k_w[1])?,
        ];
        let resistance = finite(resistances[0] + resistances[1])?;
        let bulk = finite(resistance.recip())?;
        if bulk_resistances.iter().any(|r| !positive(*r)) || !positive(bulk) {
            return Err(DepositionError::Arithmetic);
        }
        let weights = [resistances[0] / resistance, resistances[1] / resistance];
        let delta = finite(temperatures[0] - temperatures[1])?;
        // x=(Tsa-Tsb)/(Ta-Tb); never subtract rounded surface temperatures or
        // divide by a small bulk delta. Positive passivity brackets x in [0,1].
        let mut x = fraction(bulk, self.coefficients(temperatures)?.0);
        let (mut lo, mut hi) = (0.0, 1.0);
        for iterations in 1..=self.controls.maximum_iterations {
            checkpoint()?;
            let surfaces = [
                (-weights[0] * (1.0 - x)).mul_add(delta, temperatures[0]),
                (weights[1] * (1.0 - x)).mul_add(delta, temperatures[1]),
            ];
            let (surface_g, slopes, radiation_g) = self.coefficients(surfaces)?;
            let g = finite(area * series(bulk, surface_g))?;
            let power = finite(g * delta)?;
            let surface_delta = x * delta;
            let gap = finite(area * self.gap_htc_w_m2_k * surface_delta)?;
            let radiation = finite(area * radiation_g * surface_delta)?;
            let surface_power = finite(gap + radiation)?;
            let conducted = finite(area * bulk * ((1.0 - x) * delta))?;
            let residual = finite(power - surface_power)?;
            let allowance = finite(
                self.controls.absolute_power_tolerance_w
                    + self.controls.relative_power_tolerance
                        * power.abs().max(surface_power.abs()).max(conducted.abs()),
            )?;
            let scale = bulk.max(slopes[0]).max(slopes[1]).max(surface_g);
            let denominator =
                bulk / scale + weights[0] * (slopes[0] / scale) + weights[1] * (slopes[1] / scale);
            if residual.abs().max((conducted - surface_power).abs()) <= allowance {
                // q=f(Ta-ra*q,Tb+rb*q); ra'=-(da/ka)*(ka'/ka).
                // Differentiate this original balance, including the two k(h)
                // slopes even when either caloric temperature slope is zero.
                let full_drop = fraction(surface_g, bulk) * delta;
                let conductivity_terms = [
                    finite(
                        (bulk_resistances[0] / resistance)
                            * full_drop
                            * (states[0].conductivity_slope_kg_w_m_k_j / conductivities[0]),
                    )?,
                    finite(
                        (bulk_resistances[1] / resistance)
                            * full_drop
                            * (states[1].conductivity_slope_kg_w_m_k_j / conductivities[1]),
                    )?,
                ];
                let derivatives = [
                    finite(
                        area * (bulk * (slopes[0] / scale) / denominator)
                            * finite(states[0].temperature_slope_kg_k_j + conductivity_terms[0])?,
                    )?,
                    finite(
                        -area
                            * (bulk * (slopes[1] / scale) / denominator)
                            * finite(states[1].temperature_slope_kg_k_j - conductivity_terms[1])?,
                    )?,
                ];
                let link = ImplicitLinkState {
                    conductance_w_k: g,
                    power_enthalpy_derivatives_w_kg_j: derivatives,
                };
                link.admit()?;
                checkpoint()?;
                return Ok(GrayRadiativeResponse {
                    link,
                    surface_temperatures_k: Some(surfaces),
                    power_w: power,
                    gap_power_w: gap,
                    radiation_power_w: radiation,
                    surface_residual_w: residual,
                    iterations,
                });
            }
            let defect = (bulk / scale) * (1.0 - x) - (surface_g / scale) * x;
            if defect > 0.0 {
                lo = x;
            } else {
                hi = x;
            }
            let proposed = x + defect / denominator;
            let next = if proposed > lo && proposed < hi && proposed != x {
                proposed
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == x || next <= lo || next >= hi {
                break;
            }
            x = next;
        }
        Err(DepositionError::InvalidInput(
            "gray link surface balance did not converge",
        ))
    }

    fn state(
        &self,
        side: usize,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitConductionState, DepositionError> {
        checkpoint()?;
        let [lo, hi] = self.materials[side].enthalpy_bounds_j_kg();
        if !lo.is_finite() || !hi.is_finite() || !h.is_finite() || !(lo..=hi).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        let state = self.materials[side].state_at_specific_enthalpy(h)?;
        if !positive(state.temperature_k)
            || !nonnegative(state.temperature_slope_kg_k_j)
            || !nonnegative(state.conductivity_w_m_k)
            || !state.conductivity_slope_kg_w_m_k_j.is_finite()
            || (state.conductivity_w_m_k == 0.0 && state.conductivity_slope_kg_w_m_k_j != 0.0)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(state)
    }

    fn coefficients(&self, t: [f64; 2]) -> Result<(f64, [f64; 2], f64), DepositionError> {
        let mut radiation = 0.0;
        let mut slopes = [self.gap_htc_w_m2_k; 2];
        if self.effective_emissivity != 0.0 {
            let scale = self.effective_emissivity * SIGMA;
            radiation = finite(scale * (t[0] + t[1]) * t[0].mul_add(t[0], t[1] * t[1]))?;
            for side in 0..2 {
                slopes[side] =
                    finite(slopes[side] + ((4.0 * scale * t[side]) * t[side]) * t[side])?;
            }
        }
        Ok((finite(self.gap_htc_w_m2_k + radiation)?, slopes, radiation))
    }
}

impl ImplicitLinkLaw for GrayRadiativeLink<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [
            self.materials[0].caloric_law(),
            self.materials[1].caloric_law(),
        ]
    }
    fn state_at_specific_enthalpies(
        &self,
        enthalpies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        self.evaluate(enthalpies_j_kg, geometry, patch, checkpoint)
            .map(|r| r.link)
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn finite(x: f64) -> Result<f64, DepositionError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
// Positive series and fractional coefficients without overflowing products/sums.
fn series(a: f64, b: f64) -> f64 {
    if a <= b {
        a / (1.0 + a / b)
    } else {
        b / (1.0 + b / a)
    }
}
fn fraction(a: f64, b: f64) -> f64 {
    if a <= b {
        let r = a / b;
        r / (1.0 + r)
    } else {
        1.0 / (1.0 + b / a)
    }
}
