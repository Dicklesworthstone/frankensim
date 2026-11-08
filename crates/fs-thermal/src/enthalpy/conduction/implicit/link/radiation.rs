//! Two finite bodies exchanging heat through one pair of gray gap surfaces.

use super::super::boundary::{BoundGraySurface, ImplicitGraySurfaceLaw};
use super::super::contact::{ImplicitContactLaw, ImplicitContactState};
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
/// `with_surface_laws` replaces fixed emissivities with the original ordered
/// bulk-enthalpy surface sources, including their temperature/phase response.
/// `with_gap_resistance` can replace the gap HTC with an original two-sided
/// nonradiative contact law at the ordered bulk enthalpies.
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
/// enclosure solver, profile-integrated k(T), surface-temperature-dependent emissivity,
/// near-field radiation or participating-medium transport is introduced.
/// The face kernel allocates nothing; a supplied source law retains its own
/// evaluation costs and is queried once per enthalpy trial. The geometry-aware
/// trait query includes each cell's density-driven half-thickness in the same
/// balance; area, view, coatings and independently authored gap remain fixed.
pub struct GrayRadiativeLink<'a> {
    materials: [&'a dyn ImplicitConductionLaw; 2],
    emissivities: [f64; 2],
    effective_emissivity: f64,
    surface_laws: Option<[BoundGraySurface<'a>; 2]>,
    gap_htc_w_m2_k: f64,
    gap_resistance: Option<BoundGapResistance<'a>>,
    coating_resistances_m2_k_w: [f64; 2],
    controls: GrayRadiativeControls,
}

struct BoundGapResistance<'a> {
    law: &'a dyn ImplicitContactLaw,
    bounds_j_kg: [[f64; 2]; 2],
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
            emissivities,
            effective_emissivity,
            surface_laws: None,
            gap_htc_w_m2_k,
            gap_resistance: None,
            coating_resistances_m2_k_w,
            controls,
        })
    }

    /// Replace both zero fixed emissivities with original enthalpy surface laws.
    ///
    /// Each source uses the SAME ordered caloric object as its finite body.
    /// Its emissivity and complete d(epsilon)/dh are resolved once per trial
    /// at that body's original BULK enthalpy, before any insulation shortcut.
    /// Solved face temperatures do not replace the source coordinates. This is
    /// the same declared cell-state gray closure as the ambient boundary, not
    /// a model for separate surface temperatures or evolving surface phases.
    ///
    /// Both source convection values and derivatives must be zero: the common
    /// gap HTC or `with_gap_resistance` already owns nonradiative exchange.
    /// Fixed gap/coating inputs remain unchanged. Attachment is unique and
    /// requires both original fixed emissivities to be zero. Source values,
    /// domains and tangents are checked at every actual native trial.
    ///
    /// The full effective-emissivity derivative includes an opening face at
    /// epsilon=0 without dividing by that emissivity. When BOTH faces have
    /// zero emissivity and nonzero onset derivatives, the radiosity reduction
    /// has no unique joint linearization; that source corner refuses rather
    /// than fabricating a Newton tangent.
    pub fn with_surface_laws(
        mut self,
        laws: [&'a dyn ImplicitGraySurfaceLaw; 2],
    ) -> Result<Self, DepositionError> {
        if self.surface_laws.is_some() || self.emissivities != [0.0; 2] {
            return Err(DepositionError::InvalidInput(
                "gray link surface sources must uniquely replace two zero fixed emissivities",
            ));
        }
        self.surface_laws = Some([
            BoundGraySurface::new(self.materials[0].caloric_law(), laws[0])?,
            BoundGraySurface::new(self.materials[1].caloric_law(), laws[1])?,
        ]);
        Ok(self)
    }

    /// Compose one original NONRADIATIVE gap resistance with gray radiation.
    ///
    /// The source's area-specific R(h_A,h_B) is evaluated at the actual ordered
    /// BULK enthalpies, using its own temperature/phase chain rules. The solved
    /// surface temperatures do not replace those source coordinates. Pressure,
    /// separation and other external source conditions remain explicitly fixed
    /// by that law. Its entire resistance must describe only the conducting
    /// gap: do not pass a total transfer coefficient that already includes
    /// radiation, bulk half cells or the explicit coating resistances.
    ///
    /// The constructor's gap HTC must be zero and attachment may occur once.
    /// Both caloric objects must match the same ordered original materials.
    /// Actual source states are probed by native link attachment and every
    /// subsequent trial. Each must have strictly positive finite R; perfect
    /// contact, extrapolation and unrepresentable reciprocals refuse explicitly.
    /// Attach the composed link only through `with_link_laws`, not separately
    /// through `with_contacts`. Both full dR/dh contributions enter the original
    /// two-surface tangent, including latent-phase resistance feedback.
    pub fn with_gap_resistance(
        mut self,
        law: &'a dyn ImplicitContactLaw,
    ) -> Result<Self, DepositionError> {
        if self.gap_resistance.is_some() || self.gap_htc_w_m2_k != 0.0 {
            return Err(DepositionError::InvalidInput(
                "gray gap source must uniquely replace zero authored gap HTC",
            ));
        }
        let caloric = law.caloric_laws();
        let bounds_j_kg = law.enthalpy_bounds_j_kg();
        for side in 0..2 {
            let [lo, hi] = bounds_j_kg[side];
            if !lo.is_finite()
                || !hi.is_finite()
                || lo >= hi
                || !std::ptr::addr_eq(
                    self.materials[side].caloric_law() as *const dyn ThermalLaw,
                    caloric[side] as *const dyn ThermalLaw,
                )
            {
                return Err(DepositionError::InvalidInput(
                    "gray gap source requires finite bounds and original ordered caloric objects",
                ));
            }
        }
        self.gap_resistance = Some(BoundGapResistance { law, bounds_j_kg });
        Ok(self)
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
        self.evaluate_with_geometry_tangent(enthalpies_j_kg, geometry, patch, [0.0; 2], checkpoint)
    }

    fn evaluate_with_geometry_tangent(
        &self,
        enthalpies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        relative_distance_slopes_kg_j: [f64; 2],
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
            || relative_distance_slopes_kg_j
                .iter()
                .any(|slope| !slope.is_finite())
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
        // Query each original source once, even if a bulk proves insulating.
        // Neither shortcuts nor surface Newton may hide source-domain failures
        // or replace the declared bulk-state coordinates.
        let (emissivity, emissivity_slopes) = self.emissivity_at(enthalpies_j_kg, checkpoint)?;
        let (gap_htc, gap_state) = self.gap_at(enthalpies_j_kg, checkpoint)?;
        let disconnected = emissivity == 0.0 && gap_htc == 0.0;
        if (disconnected && emissivity_slopes == [0.0; 2]) || conductivities.contains(&0.0) {
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
        let mut x = fraction(
            bulk,
            self.coefficients(temperatures, gap_htc, emissivity)?.0,
        );
        let (mut lo, mut hi) = (0.0, 1.0);
        for iterations in 1..=self.controls.maximum_iterations {
            checkpoint()?;
            let surfaces = [
                (-weights[0] * (1.0 - x)).mul_add(delta, temperatures[0]),
                (weights[1] * (1.0 - x)).mul_add(delta, temperatures[1]),
            ];
            let (surface_g, slopes, radiation_g) =
                self.coefficients(surfaces, gap_htc, emissivity)?;
            let g = finite(area * series(bulk, surface_g))?;
            let power = finite(g * delta)?;
            let surface_delta = x * delta;
            let gap = finite(area * gap_htc * surface_delta)?;
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
                // q=f(Ta-ra*q,Tb+rb*q); r_i'=(d_i/k_i)*(d_i'/d_i-k_i'/k_i).
                // Differentiate this same balance, including original density
                // and conductivity changes even on an isothermal plateau.
                let full_drop = fraction(surface_g, bulk) * delta;
                let resistance_terms = [
                    finite(
                        (bulk_resistances[0] / resistance)
                            * full_drop
                            * (states[0].conductivity_slope_kg_w_m_k_j / conductivities[0]
                                - relative_distance_slopes_kg_j[0]),
                    )?,
                    finite(
                        (bulk_resistances[1] / resistance)
                            * full_drop
                            * (states[1].conductivity_slope_kg_w_m_k_j / conductivities[1]
                                - relative_distance_slopes_kg_j[1]),
                    )?,
                ];
                let mut derivatives = [
                    finite(
                        area * (bulk * (slopes[0] / scale) / denominator)
                            * finite(states[0].temperature_slope_kg_k_j + resistance_terms[0])?,
                    )?,
                    finite(
                        -area
                            * (bulk * (slopes[1] / scale) / denominator)
                            * finite(states[1].temperature_slope_kg_k_j - resistance_terms[1])?,
                    )?,
                ];
                if let Some(state) = gap_state {
                    // At fixed surfaces dQ/dR = -Q_gap/R. Their response gives
                    // D=1+ra*F'_a+rb*F'_b, so the additional full derivative is
                    // -Q_gap * R'_i / (R*D). This normalized 1/(R*D) avoids
                    // forming R^2, a reciprocal derivative, or an overflowing D.
                    let factor = finite(bulk * (gap_htc / scale) / denominator)?;
                    for side in 0..2 {
                        let feedback = finite(
                            gap * finite(
                                factor * state.resistance_enthalpy_derivatives_m2_kg_k_w_j[side],
                            )?,
                        )?;
                        derivatives[side] = finite(derivatives[side] - feedback)?;
                    }
                }
                // At fixed faces, dq/dh = sigma*(Tsa^4-Tsb^4)*epsilon_eff'.
                // The same implicit surface response divides this by
                // D=1+ra*F'_a+rb*F'_b. Factoring the fourth-power difference
                // preserves near equilibrium and the zero-emissivity onset.
                let inverse_surface_response = (bulk / scale) / denominator;
                for side in 0..2 {
                    if emissivity_slopes[side] != 0.0 {
                        let coefficient = finite(
                            (SIGMA * emissivity_slopes[side])
                                * (surfaces[0] + surfaces[1])
                                * surfaces[0].mul_add(surfaces[0], surfaces[1] * surfaces[1]),
                        )?;
                        let feedback = finite(
                            (area * inverse_surface_response)
                                * finite(coefficient * surface_delta)?,
                        )?;
                        derivatives[side] = finite(derivatives[side] + feedback)?;
                    }
                }
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

    fn gap_at(
        &self,
        enthalpies: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(f64, Option<ImplicitContactState>), DepositionError> {
        let Some(bound) = &self.gap_resistance else {
            return Ok((self.gap_htc_w_m2_k, None));
        };
        for (h, [lo, hi]) in enthalpies.into_iter().zip(bound.bounds_j_kg) {
            if !h.is_finite() || !(lo..=hi).contains(&h) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        checkpoint()?;
        let state = bound.law.state_at_specific_enthalpies(enthalpies)?;
        state.admit()?;
        if !positive(state.resistance_m2_k_w) {
            return Err(DepositionError::MaterialDomain);
        }
        let htc = finite(state.resistance_m2_k_w.recip())?;
        if !positive(htc) {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok((htc, Some(state)))
    }

    fn emissivity_at(
        &self,
        enthalpies: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(f64, [f64; 2]), DepositionError> {
        let Some(laws) = &self.surface_laws else {
            return Ok((self.effective_emissivity, [0.0; 2]));
        };
        let states = [
            laws[0].evaluate(enthalpies[0], checkpoint)?,
            laws[1].evaluate(enthalpies[1], checkpoint)?,
        ];
        if states.iter().any(|s| {
            s.convection_w_m2_k != 0.0 || s.convection_enthalpy_derivative_kg_w_m2_k_j != 0.0
        }) {
            return Err(DepositionError::MaterialDomain);
        }
        let values = states.map(|s| s.emissivity);
        let slopes = states.map(|s| s.emissivity_enthalpy_derivative_kg_j);
        let small = usize::from(values[0] > values[1]);
        let large = 1 - small;
        let low = values[small];
        let high = values[large];
        if high == 0.0 {
            if slopes.iter().all(|s| *s != 0.0) {
                return Err(DepositionError::InvalidInput(
                    "two zero-emissivity onsets have no joint gray-link tangent",
                ));
            }
            return Ok((0.0, [0.0; 2]));
        }
        // epsilon_eff=low/(1-low+low/high), d_eff/d_low=(1/denominator)^2,
        // d_eff/d_high=(low/high/denominator)^2. Every ratio is <=1; no tiny
        // emissivity reciprocal, product of emissivities or 0/0 is formed.
        let ratio = low / high;
        let denominator = (1.0 - low) + ratio;
        let low_weight = denominator.recip();
        let high_weight = ratio / denominator;
        let mut derivatives = [0.0; 2];
        derivatives[small] = finite((slopes[small] * low_weight) * low_weight)?;
        derivatives[large] = finite((slopes[large] * high_weight) * high_weight)?;
        Ok((low / denominator, derivatives))
    }

    fn coefficients(
        &self,
        t: [f64; 2],
        gap_htc: f64,
        effective_emissivity: f64,
    ) -> Result<(f64, [f64; 2], f64), DepositionError> {
        let mut radiation = 0.0;
        let mut slopes = [gap_htc; 2];
        if effective_emissivity != 0.0 {
            let scale = effective_emissivity * SIGMA;
            radiation = finite(scale * (t[0] + t[1]) * t[0].mul_add(t[0], t[1] * t[1]))?;
            for side in 0..2 {
                slopes[side] =
                    finite(slopes[side] + ((4.0 * scale * t[side]) * t[side]) * t[side])?;
            }
        }
        Ok((finite(gap_htc + radiation)?, slopes, radiation))
    }
}

impl ImplicitLinkLaw for GrayRadiativeLink<'_> {
    fn supports_current_geometry(&self) -> bool {
        true
    }

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

    fn state_at_specific_enthalpies_with_geometry_tangent(
        &self,
        enthalpies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        half_distance_relative_enthalpy_derivatives_kg_j: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitLinkState>, DepositionError> {
        self.evaluate_with_geometry_tangent(
            enthalpies_j_kg,
            geometry,
            patch,
            half_distance_relative_enthalpy_derivatives_kg_j,
            checkpoint,
        )
        .map(|response| Some(response.link))
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
