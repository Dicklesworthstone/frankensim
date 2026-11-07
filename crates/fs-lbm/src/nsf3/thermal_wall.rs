//! Spatial thermal data on the original impermeable FV wall faces.

use super::{NsfError, NsfGrid3};
use crate::lattice3::Lattice3;

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), NsfError> {
    if cancelled() {
        Err(NsfError::Input("thermal wall update cancelled"))
    } else {
        Ok(())
    }
}

/// Thermal law of one impermeable wall face. All supplied quantities are SI.
/// These are boundary conditions on the original gas energy equation, not a
/// solid-storage model or a wall-function/turbulent heat-transfer closure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NsfWallThermal {
    /// Insulated face.
    Adiabatic,
    /// Prescribed interface temperature [K], through the molecular half cell.
    Isothermal(f64),
    /// Prescribed heat flux [W/m²], positive into the gas on either side.
    HeatFlux(f64),
    /// Molecular half-cell conduction in series with a prescribed ambient film.
    /// Conductivity is evaluated at the fluid/interface mean temperature and
    /// adjacent gas pressure, using the same source as an isothermal wall.
    Film {
        coefficient_w_m2_k: f64,
        ambient_temperature_k: f64,
    },
}

/// One cell-centered gas control volume's impermeable face. `high` selects the
/// positive direction of `axis`; a neighboring solid or a non-periodic domain
/// edge must own that face. Different faces of the same cell are independent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfThermalWallFace {
    pub cell: [usize; 3],
    pub axis: usize,
    pub high: bool,
    pub thermal: NsfWallThermal,
}

impl NsfThermalWallFace {
    fn key(&self) -> ([usize; 3], usize, bool) {
        (self.cell, self.axis, self.high)
    }

    pub(super) fn feed_hash(&self, feed: &mut impl FnMut(f64)) {
        feed(10.0);
        for coordinate in self.cell {
            feed(coordinate as f64);
        }
        feed(self.axis as f64);
        feed(f64::from(u8::from(self.high)));
        match self.thermal {
            NsfWallThermal::Adiabatic => feed(0.0),
            NsfWallThermal::Isothermal(t) => {
                feed(1.0);
                feed(t);
            }
            NsfWallThermal::HeatFlux(q) => {
                feed(2.0);
                feed(q);
            }
            NsfWallThermal::Film {
                coefficient_w_m2_k,
                ambient_temperature_k,
            } => {
                feed(3.0);
                feed(coefficient_w_m2_k);
                feed(ambient_temperature_k);
            }
        }
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Attach local wall thermal data; all other faces keep `NsfSpec::wall`.
    /// The same transactional admission as [`Self::set_thermal_walls`] applies.
    pub fn with_thermal_walls(mut self, walls: Vec<NsfThermalWallFace>) -> Result<Self, NsfError> {
        self.set_thermal_walls(walls)?;
        Ok(self)
    }

    /// Replace all local wall thermal data at the accepted time. An empty list
    /// restores the global default. Populations, physical state, time and ledgers
    /// are unchanged; the original current FV heat-flux history is refreshed so
    /// the next trapezoid uses the new boundary data at both time levels.
    ///
    /// Faces must be unique impermeable fluid faces; periodic fluid neighbors,
    /// reservoir/characteristic openings and solid donor cells are refused.
    /// Temperatures (including film ambient) must lie in the declared gas window;
    /// source material failures and nonfinite heat fluxes are never extrapolated.
    /// Channel dynamic SGS retains exact uniform scalar wall data where valid;
    /// heterogeneous or flux/film planes use its existing fluid-only derivative
    /// at both physical and Favre-filtered levels. Its fit is refreshed in this
    /// same transaction; native molecular faces remain the sole wall heat owner.
    pub fn set_thermal_walls(&mut self, walls: Vec<NsfThermalWallFace>) -> Result<(), NsfError> {
        self.set_thermal_walls_cancellable(walls, || false)
    }

    /// Cancellable form of [`Self::set_thermal_walls`]. Cancellation is checked
    /// before admission, after admission, and after the complete native flux
    /// evaluation. A refusal or cancellation publishes neither controls nor history.
    pub fn set_thermal_walls_cancellable(
        &mut self,
        walls: Vec<NsfThermalWallFace>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), NsfError> {
        self.set_thermal_walls_impl(walls, cancelled, false).map(|_| ())
    }

    // The original history refresh can retain its exact sparse molecular
    // face fluxes for a subsequent native conjugate trial. Public setters keep
    // their existing behavior and discard this optional interval input.
    pub(super) fn set_thermal_walls_recorded(
        &mut self,
        walls: Vec<NsfThermalWallFace>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, NsfError> {
        self.set_thermal_walls_impl(walls, cancelled, true)
    }

    fn set_thermal_walls_impl(
        &mut self,
        mut walls: Vec<NsfThermalWallFace>,
        mut cancelled: impl FnMut() -> bool,
        record: bool,
    ) -> Result<Vec<f64>, NsfError> {
        checkpoint(&mut cancelled)?;
        if walls.len() > self.n.saturating_mul(6) {
            return Err(NsfError::Input(
                "thermal wall count exceeds the grid's face count",
            ));
        }
        let [low_k, high_k] = self.spec.temperature_window_k;
        let flux_unit = self.law.heat_flux_unit_w_m2();
        if !walls.is_empty() && !(flux_unit.is_finite() && flux_unit > 0.0) {
            return Err(NsfError::Input(
                "thermal wall heat-flux scaling must be finite and positive",
            ));
        }
        for wall in &walls {
            if wall.axis >= 3 || wall.cell.iter().zip(self.spec.dims).any(|(&x, n)| x >= n) {
                return Err(NsfError::Input("thermal wall face is outside the grid"));
            }
            let i = self.index(wall.cell);
            let direction = if wall.high { 1 } else { -1 };
            if self.solid[i]
                || self.fluid_at(wall.cell, wall.axis, direction).is_some()
                || self.open_at(i, wall.axis, direction).is_some()
                || (self.open.is_some_and(|open| open.axis == wall.axis)
                    && self.shifted(wall.cell, wall.axis, direction).is_none())
            {
                return Err(NsfError::Input(
                    "thermal data require an impermeable fluid wall face",
                ));
            }
            let temperature = match wall.thermal {
                NsfWallThermal::Adiabatic => None,
                NsfWallThermal::Isothermal(t) => Some(t),
                NsfWallThermal::HeatFlux(q) => {
                    if !q.is_finite() || !(q / flux_unit).is_finite() {
                        return Err(NsfError::Input(
                            "wall heat flux must be finite in SI and lattice units",
                        ));
                    }
                    None
                }
                NsfWallThermal::Film {
                    coefficient_w_m2_k: h,
                    ambient_temperature_k: t,
                } => {
                    if !(h.is_finite() && h >= 0.0)
                        || !(h * (3.0 * self.law.temperature_ref_k) / flux_unit).is_finite()
                    {
                        return Err(NsfError::Input(
                            "wall film coefficient must be finite and nonnegative",
                        ));
                    }
                    Some(t)
                }
            };
            if let Some(t) = temperature {
                if !(t.is_finite() && t >= low_k && t <= high_k) {
                    return Err(NsfError::Window {
                        cell: i,
                        temperature_k: t,
                        low_k,
                        high_k,
                    });
                }
            }
        }
        walls.sort_unstable_by_key(NsfThermalWallFace::key);
        if walls.windows(2).any(|pair| pair[0].key() == pair[1].key()) {
            return Err(NsfError::Input("thermal wall faces must be unique"));
        }
        checkpoint(&mut cancelled)?;
        if walls == self.thermal_walls && !record {
            return Ok(Vec::new());
        }
        let old = core::mem::replace(&mut self.thermal_walls, walls);
        let history = (|| {
            let (rho, momentum) = self.macro_fields(&self.f);
            let velocity = self.velocity_field(&rho, &momentum, &self.force);
            // Dynamic scalar derivatives depend on the actual wall table.
            // Stage their new fit before refreshing this accepted state's FV
            // history; neither its transport nor its controls publish on error.
            // Momentum, material and open-face traces stay with their owners.
            let transport = if self.dynamic_sgs.is_some() {
                Some(self.sgs_field(
                    &rho,
                    &velocity,
                    &self.theta,
                    &self.material,
                    self.open_faces.as_ref(),
                )?)
            } else {
                None
            };
            let coefficients = transport.as_deref().unwrap_or(&self.sgs_transport);
            let stress = self.stress_field(
                &self.f,
                &self.force,
                &rho,
                &self.theta,
                &self.material,
                coefficients,
            );
            let history = self.flux_divergence_with_wall_heat(
                &self.e,
                &rho,
                &velocity,
                &self.theta,
                &stress,
                coefficients,
                self.open_faces.as_ref(),
            )?;
            checkpoint(&mut cancelled)?;
            Ok((history, transport))
        })();
        let ((div, wall, open, local_heat), transport) = match history {
            Ok(history) => history,
            Err(error) => {
                self.thermal_walls = old;
                return Err(error);
            }
        };
        if let Some(transport) = transport {
            self.sgs_transport = transport;
        }
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(local_heat)
    }

    /// Current sparse overrides in canonical `(cell, axis, high)` order.
    #[must_use]
    pub fn thermal_walls(&self) -> &[NsfThermalWallFace] {
        &self.thermal_walls
    }

    pub(super) fn local_thermal_wall(
        &self,
        i: usize,
        axis: usize,
        direction: i64,
    ) -> Option<NsfWallThermal> {
        self.local_thermal_wall_index(i, axis, direction)
            .map(|j| self.thermal_walls[j].thermal)
    }

    pub(super) fn local_thermal_wall_index(
        &self,
        i: usize,
        axis: usize,
        direction: i64,
    ) -> Option<usize> {
        if self.thermal_walls.is_empty() {
            return None;
        }
        let key = (self.at(i), axis, direction > 0);
        self.thermal_walls
            .binary_search_by_key(&key, NsfThermalWallFace::key)
            .ok()
    }

    pub(super) fn thermal_wall_fluxes(
        &self,
        rho: &[f64],
        theta: &[f64],
    ) -> Result<Vec<f64>, NsfError> {
        // One sparse surface pass, no full-cell sweep. Each nonlinear film solve
        // is shared exactly by the original divergence and boundary ledger.
        let mut flux = vec![Ok(0.0); self.thermal_walls.len()];
        super::par_fill(&mut flux, self.spec.workers, &|j| {
            let wall = self.thermal_walls[j];
            let i = self.index(wall.cell);
            self.thermal_wall_flux(
                i,
                wall.axis,
                if wall.high { 1 } else { -1 },
                wall.thermal,
                rho[i],
                theta[i],
            )
        });
        flux.into_iter().collect()
    }

    pub(super) fn thermal_walls_on_domain_face(&self, axis: usize, high: bool) -> bool {
        let edge = if high { self.spec.dims[axis] - 1 } else { 0 };
        self.thermal_walls
            .iter()
            .any(|wall| wall.axis == axis && wall.high == high && wall.cell[axis] == edge)
    }

    /// Original FV sign convention: heat flux in the positive coordinate direction.
    pub(super) fn thermal_wall_flux(
        &self,
        i: usize,
        axis: usize,
        direction: i64,
        thermal: NsfWallThermal,
        rho: f64,
        theta: f64,
    ) -> Result<f64, NsfError> {
        let inward = match thermal {
            NsfWallThermal::Adiabatic => return Ok(0.0),
            NsfWallThermal::Isothermal(t) => {
                let tw = self.law.theta_of_kelvin(t);
                let k = self.spec.gas.face_conductivity(
                    &self.law,
                    0.5 * (theta + tw),
                    self.law.pressure_pa(rho, theta),
                )?;
                // Keep the original isothermal wall arithmetic on both sides.
                let flux = if direction > 0 {
                    -k * (tw - theta) / 0.5
                } else {
                    -k * (theta - tw) / 0.5
                };
                return if flux.is_finite() {
                    Ok(flux)
                } else {
                    Err(NsfError::NonFinite { cell: i })
                };
            }
            NsfWallThermal::HeatFlux(q) => q / self.law.heat_flux_unit_w_m2(),
            NsfWallThermal::Film {
                coefficient_w_m2_k,
                ambient_temperature_k,
            } => self.film_heat_flux(
                i,
                axis,
                direction > 0,
                rho,
                theta,
                coefficient_w_m2_k,
                ambient_temperature_k,
            )?,
        };
        if !inward.is_finite() {
            return Err(NsfError::NonFinite { cell: i });
        }
        Ok(if inward == 0.0 {
            0.0
        } else {
            -(direction as f64) * inward
        })
    }

    fn film_heat_flux(
        &self,
        i: usize,
        axis: usize,
        high: bool,
        rho: f64,
        theta: f64,
        coefficient: f64,
        ambient: f64,
    ) -> Result<f64, NsfError> {
        let difference = self.law.theta_of_kelvin(ambient) - theta;
        if coefficient == 0.0 || difference == 0.0 {
            return Ok(0.0);
        }
        let film =
            coefficient * (3.0 * self.law.temperature_ref_k) / self.law.heat_flux_unit_w_m2();
        let pressure = self.law.pressure_pa(rho, theta);
        let evaluate = |fraction: f64| -> Result<(f64, f64), NsfError> {
            // Tw = Tcell + fraction*(Tambient-Tcell). The source is evaluated
            // at the actual fluid/interface mean, not at bulk/ambient mean.
            let k = self.spec.gas.face_conductivity(
                &self.law,
                theta + 0.5 * fraction * difference,
                pressure,
            )?;
            let molecular = 2.0 * k;
            if !(molecular.is_finite() && molecular > 0.0) {
                return Err(NsfError::NonFinite { cell: i });
            }
            // Stable h/(2k+h) and (2k*h)/(2k+h), including almost-insulating
            // and almost-isothermal films. No large coefficient products or
            // subtraction of nearly identical interface temperatures.
            let (target, conductance) = if molecular <= film {
                let ratio = molecular / film;
                (1.0 / (1.0 + ratio), molecular / (1.0 + ratio))
            } else {
                let ratio = film / molecular;
                (ratio / (1.0 + ratio), film / (1.0 + ratio))
            };
            Ok((fraction - target, conductance * difference))
        };
        let (mut low, mut upper) = (0.0, 1.0);
        let (mut f_low, low_heat) = evaluate(low)?;
        if f_low.abs() <= 32.0 * f64::EPSILON {
            return Ok(low_heat);
        }
        let (mut f_high, high_heat) = evaluate(upper)?;
        if f_high.abs() <= 32.0 * f64::EPSILON {
            return Ok(high_heat);
        }
        for pass in 0..192 {
            // Safeguarded secants are normally a few source evaluations.
            // Every third pass bisects, bounding the work even for stiff k(T).
            let midpoint = low + 0.5 * (upper - low);
            let secant = low + (upper - low) * (-f_low / (f_high - f_low));
            let fraction = if pass % 3 == 2 || !(secant > low && secant < upper) {
                midpoint
            } else {
                secant
            };
            let (residual, inward) = evaluate(fraction)?;
            if residual.abs() <= 32.0 * f64::EPSILON {
                return Ok(inward);
            }
            if residual < 0.0 {
                low = fraction;
                f_low = residual;
            } else {
                upper = fraction;
                f_high = residual;
            }
        }
        Err(NsfError::ThermalWallConvergence {
            cell: i,
            axis,
            high,
        })
    }
}

#[cfg(test)]
mod tests;
