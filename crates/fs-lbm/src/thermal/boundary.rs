//! Thermal laws on the existing stationary, impermeable halfway slab walls.
//!
//! The scalar may be absolute temperature or a declared affine temperature
//! coordinate. Heat flux is an inward scalar increment per wall link and time
//! step. For `theta = (T - T0) / delta_T`, dimensional inputs scale as
//! `q = q'' dt / (rho Cp dx delta_T)` and `transfer = h dt / (rho Cp dx)`.
//! The original collision owns diffusivity; no relaxation factor is applied
//! again to the incoming-minus-outgoing boundary exchange.

/// Prepared spatial solid conduction coupled to the original wall reflections.
pub mod conduction;

/// One thermal law on a planar slab wall or one of its surface links.
///
/// All temperatures use the same coordinates as the thermal populations.
/// Positive heat flux enters the fluid at either wall. Geometry, material
/// properties, source scaling and validity remain with the existing slab owner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalWall {
    /// Prescribed halfway wall temperature; original anti-bounce-back law.
    Temperature(f64),
    /// Exactly zero normal scalar exchange, including nonequilibrium states.
    Adiabatic,
    /// Prescribed inward scalar increment per link over one lattice step.
    HeatFlux(f64),
    /// A film to a prescribed ambient, `q = transfer * (ambient - wall)`.
    /// The actual halfway wall temperature is solved with the reflected pair.
    Convective {
        /// Finite nonnegative `h dt / (rho Cp dx)`; zero is insulation.
        transfer: f64,
        /// Ambient temperature in the populations' scalar coordinates.
        ambient: f64,
    },
    /// A finite solid cell behind the halfway interface. Each normal link
    /// owns a separate cell; backward-Euler storage and wall conduction use
    /// that link's actual reflected heat.
    FiniteSolid(FiniteThermalWall),
}

/// One stationary solid wall cell coupled to one fluid boundary link.
///
/// Temperature uses the same affine coordinate as the thermal populations.
/// With cell heat-capacity scale `C_f = rho_f Cp_f dx^D`, the dimensionless
/// capacity is `C_s/C_f`, resistance is `R_s C_f/dt`, and prescribed heat is
/// `Q_s/(C_f delta_T)` per lattice step. `R_s` [K/W] connects the solid bulk
/// temperature to its fluid-facing surface, including declared contact
/// resistance but excluding the fluid's half cell. Zero resistance is
/// perfect contact with this lumped cell; infinite sentinels are not accepted.
///
/// Adjacent wall cells exchange no lateral heat unless a prepared wall
/// conduction graph is attached to the slab. The outside is insulated
/// except for `heat_per_step`, held over each original fluid step. This is
/// first-order lumped storage, not a resolved solid diffusion solver.
/// At extreme capacity ratios, a solid temperature increment may be below
/// floating-point resolution. Checked source/material steps reject excessive
/// independently measured first-law defects; unchecked steps report them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteThermalWall {
    /// Retained solid bulk temperature in the populations' coordinates.
    pub temperature: f64,
    /// Positive solid capacity divided by one fluid cell's heat capacity.
    pub heat_capacity: f64,
    /// Nonnegative through-wall/contact resistance in lattice units.
    pub resistance: f64,
    /// External heat into THIS solid cell per step. It is not deposited
    /// directly into the fluid a second time.
    pub heat_per_step: f64,
}

impl FiniteThermalWall {
    /// Validate the authored state before installing it on a slab.
    pub fn validate(self) -> Result<(), ThermalWallError> {
        if !self.temperature.is_finite() || !self.heat_capacity.is_finite()
            || self.heat_capacity <= 0.0 || !self.resistance.is_finite()
            || self.resistance < 0.0 || !self.heat_per_step.is_finite()
        {
            return Err(ThermalWallError::Invalid(
                "finite solid needs finite temperature/source, positive capacity and nonnegative resistance",
            ));
        }
        Ok(())
    }

    fn advance(&mut self, outgoing: f64, weight: f64) -> WallReflection {
        // q = g_in - g_out, T_s = (g_in + g_out)/(2w),
        // T_b(new) - T_s = R q, C (T_b(new) - T_b(old)) = Q - q.
        // Elimination gives q = (C (T_b(old)-g_out/w) + Q)
        // / (1 + C (R + 1/(2w))). Scaled branches avoid C*R overflow and
        // retain tiny capacities without forming Q/C in that limit.
        let fluid_trace = outgoing / weight;
        let series = self.resistance + 0.5 / weight;
        let capacity = self.heat_capacity;
        let source = self.heat_per_step;
        let (old_weight, fluid_weight, conductance, source_heat, source_temperature) =
            if capacity <= 1.0 / series {
                let ratio = capacity * series;
                let fluid_weight = 1.0 / (1.0 + ratio);
                (ratio * fluid_weight, fluid_weight, capacity * fluid_weight,
                    source * fluid_weight, (series * fluid_weight) * source)
            } else {
                let ratio = (1.0 / capacity) / series;
                let old_weight = 1.0 / (1.0 + ratio);
                (old_weight, ratio * old_weight, old_weight / series,
                    (source / capacity) * (old_weight / series),
                    (source / capacity) * old_weight)
            };
        let old = self.temperature;
        let difference = old - fluid_trace;
        let conductive_heat = if difference.is_finite() {
            conductance * difference
        } else {
            conductance * old - conductance * fluid_trace
        };
        let incoming = outgoing + (conductive_heat + source_heat);
        let heat_in = incoming - outgoing;
        // Do not update with (Q-q)/C: subtraction erases retained energy
        // when a tiny solid passes nearly all Q into the fluid. Centering
        // the convex endpoint on its dominant temperature fixes equilibrium.
        let retained = if difference.is_finite() {
            if old_weight <= fluid_weight {
                old_weight.mul_add(difference, fluid_trace)
            } else {
                (-fluid_weight).mul_add(difference, old)
            }
        } else {
            old_weight.mul_add(old, fluid_weight * fluid_trace)
        };
        self.temperature = retained + source_temperature;
        let stored_change = capacity * (self.temperature - old);
        let defect = stored_change + (heat_in - source);
        WallReflection {
            incoming,
            temperature: (0.5 * incoming + 0.5 * outgoing) / weight,
            solid: Some(FiniteThermalWallExchange {
                heat_in, stored_change, source, conducted_in: 0.0,
                maximum_link_defect: defect.abs(),
            }),
        }
    }
}

/// First-law account for finite solid links on one wall plane, in the
/// fluid's original cell-energy/temperature scale.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FiniteThermalWallExchange {
    /// Actual heat INTO the fluid from finite solid links only.
    pub heat_in: f64,
    /// Independently measured `C * (T_new - T_old)` in the solid cells.
    pub stored_change: f64,
    /// Prescribed external heat deposited into the solid cells.
    pub source: f64,
    /// Net heat into these solid cells through authored solid-conduction links.
    /// Internal graph transfers cancel across the complete solid account.
    pub conducted_in: f64,
    /// Largest absolute per-link `storage + heat_into_fluid - source - conducted_in`.
    /// Opposite local defects cannot hide in the plane sum.
    pub maximum_link_defect: f64,
}

impl FiniteThermalWallExchange {
    /// Summed first-law defect; solid/fluid exchange cancels in a combined
    /// system account, leaving only the authored external source.
    #[must_use]
    pub fn defect(self) -> f64 {
        self.stored_change + (self.heat_in - self.source) - self.conducted_in
    }
}

pub(crate) struct WallReflection {
    pub(crate) incoming: f64,
    pub(crate) temperature: f64,
    pub(crate) solid: Option<FiniteThermalWallExchange>,
}

/// Borrowed laws on one wall, in x-fastest surface order (`x` in 2-D,
/// `x + nx*y` in 3-D). Each spatial entry owns exactly one normal link.
/// Uniform fields require no per-link storage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalWallField<'a> {
    /// One law over the complete wall.
    Uniform(ThermalWall),
    /// One law per surface link; lengths are admitted by the slab owner.
    Spatial(&'a [ThermalWall]),
}

impl<'a> ThermalWallField<'a> {
    /// Authored laws, visiting a uniform declaration once. A spatial field
    /// visits every surface link, even if several links use the same law.
    pub fn laws(self) -> impl Iterator<Item = ThermalWall> + 'a {
        let (uniform, spatial) = match self {
            Self::Uniform(wall) => (Some(wall), &[][..]),
            Self::Spatial(walls) => (None, walls),
        };
        uniform.into_iter().chain(spatial.iter().copied())
    }

    /// The common law, if this field is spatially uniform. This also detects
    /// constant spatial inputs so scalar diagnostics and SGS retain their
    /// original fixed-temperature or insulating boundary treatment.
    #[must_use]
    pub fn uniform(self) -> Option<ThermalWall> {
        match self {
            Self::Uniform(wall) => Some(wall),
            Self::Spatial(walls) => walls.first().copied()
                .filter(|first| walls.iter().all(|wall| wall == first)),
        }
    }

    pub(crate) fn validate(self, links: usize) -> Result<(), ThermalWallError> {
        if let Self::Spatial(walls) = self
            && (links == 0 || walls.len() != links)
        {
            return Err(ThermalWallError::Invalid("one thermal law per wall link is required"));
        }
        for wall in self.laws() { wall.validate()?; }
        Ok(())
    }

    pub(crate) fn at(self, link: usize) -> ThermalWall {
        match self {
            Self::Uniform(wall) => wall,
            Self::Spatial(walls) => walls[link],
        }
    }

    pub(crate) fn pair(
        uniform: [ThermalWall; 2], spatial: Option<&'a [Vec<ThermalWall>; 2]>,
    ) -> [Self; 2] {
        spatial.map_or_else(
            || uniform.map(Self::Uniform),
            |walls| [Self::Spatial(&walls[0]), Self::Spatial(&walls[1])],
        )
    }
}

/// Invalid wall data; admission precedes mutation of the original slab.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermalWallError {
    /// A coefficient, temperature or prescribed flux cannot be admitted.
    Invalid(&'static str),
}

impl core::fmt::Display for ThermalWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "thermal wall: {reason}"),
        }
    }
}
impl std::error::Error for ThermalWallError {}

/// Measured exchange and surface trace for the accepted reflected populations.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ThermalWallExchange {
    /// Sum of actual incoming-minus-outgoing populations; positive into fluid.
    /// Multiply by the original physical cell heat-capacity/temperature scale
    /// to obtain energy. This is not an additional source to deposit.
    pub heat_in: f64,
    /// Arithmetic mean of halfway wall temperatures over its normal links.
    /// This is the surface trace, not the neighboring fluid-cell temperature.
    /// A prescribed temperature retains its exact authored value; the other
    /// wall laws reconstruct their trace from the accepted normal pair.
    pub mean_temperature: f64,
    /// Smallest individual wall trace; the checked owner admits local extrema.
    pub min_temperature: f64,
    /// Largest individual wall trace; a mean alone can conceal unsupported heat.
    pub max_temperature: f64,
    /// Reciprocal storage/source account when this plane has finite cells.
    /// `heat_in` above continues to include ALL boundary links.
    pub solid: Option<FiniteThermalWallExchange>,
}

impl ThermalWall {
    /// Check the authored scalar data. Negative temperature coordinates and
    /// cooling fluxes are allowed; a physical material owner admits absolute
    /// temperature and the resolved wall trace in its own validity interval.
    pub fn validate(self) -> Result<(), ThermalWallError> {
        let invalid = ThermalWallError::Invalid;
        match self {
            Self::FiniteSolid(solid) => solid.validate(),
            Self::Temperature(t) if !t.is_finite() => {
                Err(invalid("wall temperature must be finite"))
            }
            Self::HeatFlux(q) if !q.is_finite() => {
                Err(invalid("inward wall heat flux must be finite"))
            }
            Self::Convective { transfer, ambient }
                if !transfer.is_finite() || transfer < 0.0 || !ambient.is_finite() =>
            {
                Err(invalid(
                    "film transfer must be finite nonnegative and ambient finite",
                ))
            }
            _ => Ok(()),
        }
    }

    /// Reflect one post-collision outgoing population without allocating.
    ///
    /// The original ABB pair gives `T_w = (g_in + g_out)/(2w)`, while its
    /// conservative link exchange is `q = g_in - g_out`. Eliminating `T_w`
    /// from `q = H (T_ambient - T_w)` gives
    /// `g_in = g_out + 2H/(2w + H) * (w T_ambient - g_out)`.
    /// Two scaled branches retain both the insulating and prescribed-T limits
    /// without forming an overflowing `H/(2w)` or `2w + H`. This is a local
    /// planar boundary closure, not a curved-wall or unsteady-order theorem.
    /// Inputs have been admitted by the calling owner; nonrepresentable output
    /// remains visible to its existing numerical-state/temperature checks.
    pub(crate) fn reflect(self, outgoing: f64, weight: f64) -> (f64, f64) {
        let incoming = match self {
            Self::FiniteSolid(mut solid) => {
                let response = solid.advance(outgoing, weight);
                return (response.incoming, response.temperature);
            }
            // Retain the original fixed-temperature arithmetic exactly.
            Self::Temperature(temperature) => {
                return (2.0f64.mul_add(weight * temperature, -outgoing), temperature);
            }
            Self::Adiabatic => outgoing,
            Self::HeatFlux(flux) => outgoing + flux,
            Self::Convective { transfer, .. } if transfer == 0.0 => outgoing,
            Self::Convective { transfer, ambient } => {
                let target = weight * ambient;
                let pair_weight = 2.0 * weight;
                if transfer <= pair_weight {
                    let ratio = transfer / pair_weight;
                    let fraction = (2.0 * ratio) / (1.0 + ratio);
                    let difference = target - outgoing;
                    if difference.is_finite() {
                        fraction.mul_add(difference, outgoing)
                    } else {
                        // A convex reflection can remain finite even when
                        // subtracting opposite extreme populations overflows.
                        (1.0 - fraction).mul_add(outgoing, fraction * target)
                    }
                } else {
                    let ratio = pair_weight / transfer;
                    let fraction = (2.0 * ratio) / (1.0 + ratio);
                    let prescribed = 2.0f64.mul_add(target, -outgoing);
                    let difference = outgoing - target;
                    if prescribed.is_finite() && difference.is_finite() {
                        fraction.mul_add(difference, prescribed)
                    } else {
                        ((ratio - 1.0) / (1.0 + ratio))
                            .mul_add(outgoing, (2.0 / (1.0 + ratio)) * target)
                    }
                }
            }
        };
        (incoming, (0.5 * incoming + 0.5 * outgoing) / weight)
    }

    /// Commit one reflection and its solid endpoint. The existing slab
    /// source transaction owns rollback together with fluid populations.
    pub(crate) fn reflect_and_advance(&mut self, outgoing: f64, weight: f64) -> WallReflection {
        if let Self::FiniteSolid(solid) = self {
            solid.advance(outgoing, weight)
        } else {
            let (incoming, temperature) = self.reflect(outgoing, weight);
            WallReflection { incoming, temperature, solid: None }
        }
    }
}

#[cfg(test)]
mod finite_tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-13 * (1.0 + a.abs().max(b.abs())), "{a} != {b}");
    }

    #[test]
    fn finite_solid_closes_interface_conduction_and_reciprocal_storage() {
        for weight in [1.0 / 6.0, 0.125] {
            for capacity in [0.001, 0.2, 1.0, 20.0] {
                for resistance in [0.0, 0.5, 30.0] {
                    for source in [-0.03, 0.0, 0.05] {
                        for (old, fluid) in [(1.2, 0.2), (-0.5, 0.7)] {
                            let mut wall = FiniteThermalWall {
                                temperature: old, heat_capacity: capacity,
                                resistance, heat_per_step: source,
                            };
                            let result = wall.advance(weight * fluid, weight);
                            let audit = result.solid.unwrap();
                            let exact_heat = (capacity * (old - fluid) + source)
                                / (1.0 + capacity * (resistance + 0.5 / weight));
                            close(audit.heat_in, exact_heat);
                            close(wall.temperature - result.temperature, resistance * audit.heat_in);
                            close(audit.stored_change, source - audit.heat_in);
                            close(audit.defect(), 0.0);
                            assert!(audit.maximum_link_defect < 1e-13);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn finite_capacity_limits_stay_bounded_and_conductivity_changes_transfer() {
        for weight in [1.0 / 6.0, 0.125] {
            for capacity in [1e-200, 1e-10, 1.0, 1e10, 1e200] {
                let mut fast = FiniteThermalWall {
                    temperature: 1.0, heat_capacity: capacity,
                    resistance: 0.0, heat_per_step: 0.0,
                };
                let mut slow = FiniteThermalWall { resistance: 100.0, ..fast };
                let q_fast = fast.advance(0.25 * weight, weight).solid.unwrap().heat_in;
                let q_slow = slow.advance(0.25 * weight, weight).solid.unwrap().heat_in;
                assert!((0.25..=1.0).contains(&fast.temperature));
                assert!((0.25..=1.0).contains(&slow.temperature));
                assert!(q_slow <= q_fast);
                assert!(fast.temperature <= slow.temperature);
            }
            // Tiny storage passes the source without Q/C overflow, while
            // retaining the finite conducting bulk/surface temperature drop.
            let mut small = FiniteThermalWall {
                temperature: 0.5, heat_capacity: 1e-320,
                resistance: 2.0, heat_per_step: 0.3,
            };
            let result = small.advance(weight * 0.5, weight);
            close(result.solid.unwrap().heat_in, 0.3);
            close(small.temperature, 0.5 + 0.3 * (2.0 + 0.5 / weight));
        }
    }

    #[test]
    fn unrepresentable_solid_storage_remains_an_explicit_balance_defect() {
        for weight in [1.0 / 6.0, 0.125] {
            let mut wall = FiniteThermalWall {
                temperature: 1.0,
                heat_capacity: 1e300,
                resistance: 0.0,
                heat_per_step: 0.0,
            };
            let result = wall.advance(weight * 0.25, weight);
            let audit = result.solid.unwrap();
            assert_eq!(wall.temperature, 1.0);
            assert_eq!(audit.stored_change, 0.0);
            assert!(audit.heat_in > 0.1);
            assert_eq!(audit.defect(), audit.heat_in);
            assert_eq!(audit.maximum_link_defect, audit.heat_in);
        }
    }

    #[test]
    fn equilibrium_is_fixed_and_invalid_finite_walls_refuse() {
        let base = FiniteThermalWall {
            temperature: 0.75, heat_capacity: 3.0, resistance: 2.0, heat_per_step: 0.0,
        };
        for weight in [1.0 / 6.0, 0.125] {
            let mut wall = base;
            let result = wall.advance(weight * base.temperature, weight);
            assert_eq!(wall, base);
            assert_eq!(result.solid.unwrap(), FiniteThermalWallExchange::default());
        }
        assert!(base.validate().is_ok());
        for wall in [
            FiniteThermalWall { temperature: f64::NAN, ..base },
            FiniteThermalWall { heat_capacity: 0.0, ..base },
            FiniteThermalWall { heat_capacity: f64::INFINITY, ..base },
            FiniteThermalWall { resistance: -1.0, ..base },
            FiniteThermalWall { resistance: f64::INFINITY, ..base },
            FiniteThermalWall { heat_per_step: f64::NAN, ..base },
        ] { assert!(ThermalWall::FiniteSolid(wall).validate().is_err()); }
    }
}
