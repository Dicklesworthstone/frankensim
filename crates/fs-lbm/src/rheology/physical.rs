//! Physical constitutive parameters bound to the existing material scaling.
//!
//! [`MaterialScaling::bind_rheology`] uses the original material's frozen
//! reference density and its admitted dx/dt/Mach mapping. The supplied law is
//! an explicit scenario model; a transport snapshot alone does not establish
//! non-Newtonian parameters or extend the source's applicability envelope.
//!
//! The existing [`LiquidViscosityBasis`] declares the units of every viscous
//! coefficient. Dynamic input uses viscosity [Pa s], consistency [Pa s^n],
//! and yield stress [Pa]. Kinematic input uses viscosity [m²/s], consistency
//! [m² s^(n-2)], and yield stress divided by density [m²/s²]. Time and shear
//! rate parameters always use seconds and inverse seconds, respectively.

use fs_material::fluid_transport::FluidTransport;
use fs_material::liquid_temperature::LiquidViscosityBasis;

use super::Rheology;
use crate::material::MaterialScaling;

/// Constitutive parameters in SI, with dynamic/kinematic basis selected when
/// binding to a material scale. Coefficients are never inferred from a name,
/// temperature, or the reference snapshot's one Newtonian viscosity value.
#[derive(Debug, Clone, Copy)]
pub enum PhysicalRheology {
    /// Constant dynamic [Pa s] or kinematic [m²/s] viscosity.
    Newtonian {
        /// Positive viscosity in the declared basis.
        viscosity: f64,
    },
    /// Regularized power law: viscosity = K max(gamma, gamma_min)^(n-1).
    /// The SI cutoff is explicit so changing dt cannot change the material.
    PowerLaw {
        /// Positive consistency, [Pa s^n] or [m² s^(n-2)].
        consistency: f64,
        /// Positive dimensionless flow index.
        n: f64,
        /// Positive physical low-rate cutoff [1/s].
        minimum_shear_rate_per_s: f64,
    },
    /// Carreau viscosity, using the same model as [`Rheology::Carreau`].
    Carreau {
        /// Positive zero-shear viscosity, [Pa s] or [m²/s].
        viscosity_zero: f64,
        /// Positive infinite-shear viscosity in the same basis.
        viscosity_infinite: f64,
        /// Nonnegative physical relaxation time [s].
        relaxation_time_s: f64,
        /// Positive dimensionless flow index.
        n: f64,
    },
    /// Carreau–Yasuda viscosity with an independently calibrated transition
    /// exponent, using the same model as [`Rheology::CarreauYasuda`].
    CarreauYasuda {
        /// Positive zero-shear viscosity, [Pa s] or [m²/s].
        viscosity_zero: f64,
        /// Positive infinite-shear viscosity in the same basis.
        viscosity_infinite: f64,
        /// Nonnegative physical relaxation time [s].
        relaxation_time_s: f64,
        /// Positive dimensionless flow index.
        n: f64,
        /// Positive dimensionless transition exponent; unchanged by scaling.
        a: f64,
    },
    /// Regularized Herschel–Bulkley; n=1 gives Bingham–Papanastasiou.
    HerschelBulkley {
        /// Positive consistency, [Pa s^n] or [m² s^(n-2)].
        consistency: f64,
        /// Positive dimensionless flow index.
        n: f64,
        /// Nonnegative yield stress [Pa] or yield stress / density [m²/s²].
        yield_stress: f64,
        /// Nonnegative Papanastasiou regularization time [s].
        regularization_time_s: f64,
        /// Positive physical cutoff for the power-law contribution [1/s].
        minimum_shear_rate_per_s: f64,
    },
}

/// Invalid physical input or a conversion that cannot retain a nonzero term.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhysicalRheologyError {
    /// A physical parameter violates the constitutive or scaling requirements.
    Invalid(&'static str),
    /// A nonzero converted parameter overflowed or underflowed to zero.
    Unrepresentable(&'static str),
}

impl core::fmt::Display for PhysicalRheologyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(parameter) => write!(f, "invalid physical rheology: {parameter}"),
            Self::Unrepresentable(parameter) => {
                write!(f, "unrepresentable lattice rheology: {parameter}")
            }
        }
    }
}
impl std::error::Error for PhysicalRheologyError {}

/// One converted constitutive law retaining its original source/scale binding.
///
/// The material reference density is constant for this incompressible law.
/// An evolving temperature/composition model can supply another physical law
/// and replace its converted entry through [`super::RheologyField::set_law`]
/// or refresh the complete field. No new temperature or material solver runs
/// here. Conversion does not certify the lattice tau interval over all rates:
/// the existing local updater still reports numerical floor/cap events.
pub struct ScaledRheology<'s, 'm, F: FluidTransport + ?Sized> {
    scaling: &'s MaterialScaling<'m, F>,
    physical: PhysicalRheology,
    basis: LiquidViscosityBasis,
    lattice: Rheology,
}

impl<'m, F: FluidTransport + ?Sized> MaterialScaling<'m, F> {
    /// Convert caller-declared SI rheology with this source's density and
    /// already admitted physical scale. Dynamic coefficients divide by the
    /// frozen material density; kinematic coefficients do not divide again.
    ///
    /// The lattice consistency is `(K/rho) dt^(2-n)/dx²`, yield contribution
    /// is `(stress/rho) dt²/dx²`, both times divide by dt, and the physical
    /// minimum shear rate multiplies dt. The same result serves both dimensions.
    pub fn bind_rheology(
        &self,
        physical: PhysicalRheology,
        basis: LiquidViscosityBasis,
    ) -> Result<ScaledRheology<'_, 'm, F>, PhysicalRheologyError> {
        let lattice = physical.convert(
            self.cell_size_m(),
            self.timestep_s(),
            self.transport().density_kg_m3(),
            basis,
        )?;
        Ok(ScaledRheology {
            scaling: self,
            physical,
            basis,
            lattice,
        })
    }
}

impl<'s, 'm, F: FluidTransport + ?Sized> ScaledRheology<'s, 'm, F> {
    /// Original source, frozen density and common physical/lattice mapping.
    #[must_use]
    pub const fn scaling(&self) -> &'s MaterialScaling<'m, F> {
        self.scaling
    }
    /// Original SI parameters; their coefficients retain [`Self::basis`].
    #[must_use]
    pub const fn physical(&self) -> PhysicalRheology {
        self.physical
    }
    /// Whether the supplied SI viscous coefficients were dynamic or kinematic.
    #[must_use]
    pub const fn basis(&self) -> LiquidViscosityBasis {
        self.basis
    }
    /// Converted law for uniform updates or [`super::RheologyField`]. Physical
    /// power laws use zero-yield Herschel–Bulkley to preserve their SI cutoff;
    /// the legacy lattice PowerLaw variant has a fixed 1e-12 lattice cutoff.
    /// Both D2Q9 and D3Q19 evaluate the physical three-dimensional deviatoric
    /// strain invariant; a planar flow embeds its unmodeled strain as zero.
    #[must_use]
    pub const fn lattice_law(&self) -> Rheology {
        self.lattice
    }
}

impl PhysicalRheology {
    // This pure conversion is shared by every source binding and dimension.
    // Only MaterialScaling exposes it publicly, so dx/dt/Mach admission stays
    // with the existing owner instead of introducing another scaling contract.
    fn convert(
        self,
        dx: f64,
        dt: f64,
        density: f64,
        basis: LiquidViscosityBasis,
    ) -> Result<Rheology, PhysicalRheologyError> {
        positive(dx, "cell size [m]")?;
        positive(dt, "timestep [s]")?;
        positive(density, "reference density [kg/m³]")?;
        let divisor = match basis {
            LiquidViscosityBasis::Dynamic => density,
            LiquidViscosityBasis::Kinematic => 1.0,
        };
        let viscosity = |value| coefficient(value, divisor, dx, dt, 0.0, "viscosity");
        let law = match self {
            Self::Newtonian { viscosity: value } => {
                positive(value, "viscosity")?;
                Rheology::Newtonian {
                    nu: viscosity(value)?,
                }
            }
            Self::PowerLaw {
                consistency,
                n,
                minimum_shear_rate_per_s,
            } => {
                positive(consistency, "power-law consistency")?;
                positive(n, "power-law index")?;
                positive(minimum_shear_rate_per_s, "minimum shear rate [1/s]")?;
                Rheology::HerschelBulkley {
                    k: coefficient(
                        consistency,
                        divisor,
                        dx,
                        dt,
                        1.0 - n,
                        "power-law consistency",
                    )?,
                    n,
                    yield_stress_over_density: 0.0,
                    regularization_time: 0.0,
                    minimum_shear_rate: checked(
                        minimum_shear_rate_per_s * dt,
                        "minimum shear rate",
                    )?,
                }
            }
            Self::Carreau {
                viscosity_zero,
                viscosity_infinite,
                relaxation_time_s,
                n,
            } => {
                // Parameter ordering is admitted in SI before conversion;
                // rounding cannot turn an inadmissible thickening law into
                // an apparently constant one.
                Rheology::Carreau {
                    nu0: viscosity_zero,
                    nu_inf: viscosity_infinite,
                    lambda: relaxation_time_s,
                    n,
                }
                .check_parameters()
                .map_err(PhysicalRheologyError::Invalid)?;
                Rheology::Carreau {
                    nu0: viscosity(viscosity_zero)?,
                    nu_inf: viscosity(viscosity_infinite)?,
                    lambda: time(relaxation_time_s, dt)?,
                    n,
                }
            }
            Self::CarreauYasuda {
                viscosity_zero,
                viscosity_infinite,
                relaxation_time_s,
                n,
                a,
            } => {
                Rheology::CarreauYasuda {
                    nu0: viscosity_zero,
                    nu_inf: viscosity_infinite,
                    lambda: relaxation_time_s,
                    n,
                    a,
                }
                .check_parameters()
                .map_err(PhysicalRheologyError::Invalid)?;
                Rheology::CarreauYasuda {
                    nu0: viscosity(viscosity_zero)?,
                    nu_inf: viscosity(viscosity_infinite)?,
                    lambda: time(relaxation_time_s, dt)?,
                    n,
                    a,
                }
            }
            Self::HerschelBulkley {
                consistency,
                n,
                yield_stress,
                regularization_time_s,
                minimum_shear_rate_per_s,
            } => {
                positive(consistency, "Herschel-Bulkley consistency")?;
                positive(n, "Herschel-Bulkley index")?;
                nonnegative(yield_stress, "yield stress")?;
                positive(minimum_shear_rate_per_s, "minimum shear rate [1/s]")?;
                Rheology::HerschelBulkley {
                    k: coefficient(
                        consistency,
                        divisor,
                        dx,
                        dt,
                        1.0 - n,
                        "Herschel-Bulkley consistency",
                    )?,
                    n,
                    yield_stress_over_density: if yield_stress == 0.0 {
                        0.0
                    } else {
                        coefficient(yield_stress, divisor, dx, dt, 1.0, "yield stress")?
                    },
                    regularization_time: time(regularization_time_s, dt)?,
                    minimum_shear_rate: checked(
                        minimum_shear_rate_per_s * dt,
                        "minimum shear rate",
                    )?,
                }
            }
        };
        law.check_parameters()
            .map_err(PhysicalRheologyError::Invalid)?;
        Ok(law)
    }
}

// value/divisor * dt^(1+extra_time_power)/dx². Normal intermediates use the
// existing material diffusivity operation order. A strict-log fallback avoids
// losing a representable result to an intermediate power/product/division;
// it is only used when that ordinary operation tree leaves the normal range.
fn coefficient(
    value: f64,
    divisor: f64,
    dx: f64,
    dt: f64,
    extra_time_power: f64,
    name: &'static str,
) -> Result<f64, PhysicalRheologyError> {
    let kinematic = value / divisor;
    let velocity = dt / dx;
    let diffusivity = velocity / dx;
    let power = fs_math::det::pow(dt, extra_time_power);
    let scaled = kinematic * diffusivity;
    let value_lattice = scaled * power;
    if [
        kinematic,
        velocity,
        diffusivity,
        power,
        scaled,
        value_lattice,
    ]
    .iter()
    .all(|v| v.is_normal())
    {
        return Ok(value_lattice);
    }
    let logarithm = fs_math::det::ln(value) - fs_math::det::ln(divisor)
        + (1.0 + extra_time_power) * fs_math::det::ln(dt)
        - 2.0 * fs_math::det::ln(dx);
    checked(fs_math::det::exp(logarithm), name)
}

fn time(value_s: f64, dt: f64) -> Result<f64, PhysicalRheologyError> {
    nonnegative(value_s, "constitutive time [s]")?;
    if value_s == 0.0 {
        Ok(0.0)
    } else {
        checked(value_s / dt, "constitutive time")
    }
}
fn positive(value: f64, name: &'static str) -> Result<(), PhysicalRheologyError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(PhysicalRheologyError::Invalid(name))
    }
}
fn nonnegative(value: f64, name: &'static str) -> Result<(), PhysicalRheologyError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(PhysicalRheologyError::Invalid(name))
    }
}
fn checked(value: f64, name: &'static str) -> Result<f64, PhysicalRheologyError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(PhysicalRheologyError::Unrepresentable(name))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independently synthesize incoming quadratic-equilibrium populations
    // with Pi = -2 rho cs² tau D and the Guo half-source transform. The
    // antisymmetric part of grad(u) creates no viscous stress.
    fn strained_populations<const D: usize, const Q: usize>(
        directions: [[i32; D]; Q],
        weights: [f64; Q],
        gradient_per_s: [[f64; D]; D],
        dt: f64,
        tau: f64,
    ) -> [f64; Q] {
        let velocity = [0.04, -0.03, 0.0];
        let force = [0.0002, -0.0003, 0.0];
        let rho = 1.2;
        core::array::from_fn(|q| {
            let e = directions[q].map(f64::from);
            let eu: f64 = (0..D).map(|a| e[a] * velocity[a]).sum();
            let u2: f64 = (0..D).map(|a| velocity[a] * velocity[a]).sum();
            let equilibrium = weights[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
            let mut contraction = 0.0;
            for a in 0..D {
                for b in 0..D {
                    let stress = -rho
                        * crate::CS2
                        * tau
                        * dt
                        * (gradient_per_s[a][b] + gradient_per_s[b][a]);
                    contraction += (e[a] * e[b] - if a == b { crate::CS2 } else { 0.0 }) * stress;
                }
            }
            let source: f64 = (0..D)
                .map(|a| (3.0 * (e[a] - velocity[a]) + 9.0 * eu * e[a]) * force[a])
                .sum();
            equilibrium
                + weights[q] * (contraction / (2.0 * crate::CS2 * crate::CS2) - 0.5 * source)
        })
    }

    fn models() -> [PhysicalRheology; 5] {
        [
            PhysicalRheology::Newtonian { viscosity: 0.08 },
            PhysicalRheology::PowerLaw {
                consistency: 0.04,
                n: 0.6,
                minimum_shear_rate_per_s: 0.002,
            },
            PhysicalRheology::Carreau {
                viscosity_zero: 0.4,
                viscosity_infinite: 0.02,
                relaxation_time_s: 0.03,
                n: 0.4,
            },
            PhysicalRheology::HerschelBulkley {
                consistency: 0.05,
                n: 0.7,
                yield_stress: 0.2,
                regularization_time_s: 0.6,
                minimum_shear_rate_per_s: 0.003,
            },
            PhysicalRheology::CarreauYasuda {
                viscosity_zero: 0.4,
                viscosity_infinite: 0.02,
                relaxation_time_s: 0.03,
                n: 0.4,
                a: 0.7,
            },
        ]
    }

    // Independent SI evaluation, without the lattice conversion or local solve.
    fn viscosity_si(law: PhysicalRheology, rate: f64) -> f64 {
        match law {
            PhysicalRheology::Newtonian { viscosity } => viscosity,
            PhysicalRheology::PowerLaw {
                consistency,
                n,
                minimum_shear_rate_per_s,
            } => consistency * fs_math::det::pow(rate.max(minimum_shear_rate_per_s), n - 1.0),
            PhysicalRheology::Carreau {
                viscosity_zero,
                viscosity_infinite,
                relaxation_time_s,
                n,
            } => {
                viscosity_infinite
                    + (viscosity_zero - viscosity_infinite)
                        * fs_math::det::pow(
                            1.0 + (relaxation_time_s * rate).powi(2),
                            (n - 1.0) / 2.0,
                        )
            }
            PhysicalRheology::HerschelBulkley {
                consistency,
                n,
                yield_stress,
                regularization_time_s,
                minimum_shear_rate_per_s,
            } => {
                let base =
                    consistency * fs_math::det::pow(rate.max(minimum_shear_rate_per_s), n - 1.0);
                base + if rate == 0.0 {
                    yield_stress * regularization_time_s
                } else {
                    yield_stress * (-fs_math::det::expm1(-regularization_time_s * rate)) / rate
                }
            }
            PhysicalRheology::CarreauYasuda {
                viscosity_zero,
                viscosity_infinite,
                relaxation_time_s,
                n,
                a,
            } => {
                viscosity_infinite
                    + (viscosity_zero - viscosity_infinite)
                        * (1.0 + (relaxation_time_s * rate).powf(a)).powf((n - 1.0) / a)
            }
        }
    }

    #[test]
    fn physical_stress_and_local_roots_survive_refinement_in_every_law() {
        let density = 800.0;
        for model in models() {
            for refinement in [1.0, 2.0, 4.0] {
                let dx = 0.002 / refinement;
                let dt = 0.002 / (refinement * refinement);
                let law = model
                    .convert(dx, dt, density, LiquidViscosityBasis::Dynamic)
                    .unwrap();
                for rate in [0.0, 1e-4, 0.01, 1.0, 20.0, 300.0] {
                    let mu = viscosity_si(model, rate);
                    let rate_lattice = rate * dt;
                    let expected_nu = (mu / density) * dt / (dx * dx);
                    let nu = law.viscosity(rate_lattice);
                    assert!(
                        (nu / expected_nu - 1.0).abs() < 5e-13,
                        "{model:?}, rate={rate}, refinement={refinement}"
                    );
                    let physical_stress = density * (dx / dt).powi(2) * nu * rate_lattice;
                    assert!(
                        (physical_stress - mu * rate).abs() < 5e-13 * (mu * rate).abs().max(1e-12)
                    );
                    let tau = 0.5 + expected_nu / crate::CS2;
                    assert!(tau > super::super::TAU_FLOOR && tau < super::super::TAU_CAP);
                    for (numerator, fixed) in [
                        (rate_lattice * tau, 0.0),
                        (0.6 * rate_lattice * tau, 0.8 * rate_lattice),
                    ] {
                        let solved = law.prepare().solve_with_fixed_strain(numerator, fixed, 2.7);
                        assert!((solved.tau - tau).abs() < 4e-12);
                        assert!(
                            (solved.shear_rate - rate_lattice).abs()
                                < 4e-12 * rate_lattice.max(1e-20)
                        );
                        assert!(!solved.floored && !solved.capped);
                    }
                }
            }
        }
    }

    /// G0/G3: an SI law gives the same viscosity and strain in D2Q9 and
    /// its D3Q19 planar embedding, including extension, compression, shear
    /// and rigid rotation. The old full-2-D norm overestimated the rate of
    /// equal planar extension by sqrt(3), changing its apparent viscosity.
    #[test]
    fn physical_planar_rheology_agrees_between_native_dimensions() {
        use crate::d3q19::{E3, PlatesGrid3, W3, update_tau3};
        use crate::material::plan_scaling_from_material;
        use crate::rheology::update_tau;
        use fs_material::fluid_transport::FluidTransportState;

        let material = FluidTransportState::try_new(800.0, 2000.0, 0.15, 0.08).unwrap();
        let (dx, dt) = (0.002, 0.002);
        let scaling = plan_scaling_from_material(&material, dx, dt, 0.01).unwrap();
        for gradient in [
            [[4.0, 0.0], [0.0, 4.0]],
            [[-4.0, 0.0], [0.0, -4.0]],
            [[3.0, -4.0], [-4.0, -1.0]],
            [[0.0, 12.0], [0.0, 0.0]],
            [[0.0, -4.0], [4.0, 0.0]],
        ] {
            let gradient3 = [
                [gradient[0][0], gradient[0][1], 0.0],
                [gradient[1][0], gradient[1][1], 0.0],
                [0.0; 3],
            ];
            let mean = (gradient[0][0] + gradient[1][1]) / 3.0;
            let mut rate_squared = 0.0_f64;
            for a in 0..3 {
                for b in 0..3 {
                    let strain = 0.5 * (gradient3[a][b] + gradient3[b][a]);
                    let deviator = strain - if a == b { mean } else { 0.0 };
                    rate_squared += 2.0 * deviator * deviator;
                }
            }
            let rate = rate_squared.sqrt();
            for model in models() {
                let binding = scaling
                    .bind_rheology(model, LiquidViscosityBasis::Dynamic)
                    .unwrap();
                let law = binding.lattice_law();
                let nu = viscosity_si(model, rate) / 800.0 * dt / (dx * dx);
                let expected_tau = 0.5 + nu / crate::CS2;
                assert!(
                    expected_tau > super::super::TAU_FLOOR && expected_tau < super::super::TAU_CAP
                );
                let density = 1.7;
                let mut planar = crate::Grid::uniform(2, 2, density);
                planar.f.fill(strained_populations(
                    crate::E.map(|(x, y)| [x, y]),
                    crate::W,
                    gradient,
                    dt,
                    expected_tau,
                ));
                planar.fext.fill([0.0002, -0.0003]);
                let mut spatial = PlatesGrid3::plates(1, 1, 2, density);
                spatial.f.fill(strained_populations(
                    E3.map(|(x, y, z)| [x, y, z]),
                    W3,
                    gradient3,
                    dt,
                    expected_tau,
                ));
                spatial.fext.fill([0.0002, -0.0003, 0.0]);
                let planar_stats = update_tau(&mut planar, law);
                let spatial_stats = update_tau3(&mut spatial, law);
                // The law is compared at the rate the populations carry, not at the oracle's rate. A rigid
                // rotation has exactly zero strain, so what its populations carry is roundoff (measured
                // 1.55e-14 1/s), and a law whose slope at 0 is infinite (Carreau-Yasuda with a < 1, where
                // the viscosity moves with (lambda x rate)^a) turns that roundoff into a tau offset of
                // 1.08e-11 (measured, a = 0.7): a zero-rate oracle cannot be met to 4e-12 by a correct
                // update. The rate itself is still held to the oracle, and for a rotation to roundoff size.
                let tau_at = |at_rate: f64| {
                    0.5 + viscosity_si(model, at_rate) / 800.0 * dt / (dx * dx) / crate::CS2
                };
                let mut carried = [0.0; 2];
                for (index, stats) in [planar_stats, spatial_stats].into_iter().enumerate() {
                    let carried_rate = stats.max_shear / dt;
                    carried[index] = carried_rate;
                    assert!(
                        (carried_rate - rate).abs() < 2e-10,
                        "{model:?}, {gradient:?}: {carried_rate} != {rate}"
                    );
                    if rate == 0.0 {
                        // Roundoff size: at most c = 16 units of eps times the population scale (the
                        // density), 16 x 2.2e-16 x 1.7 = 6.0e-15 per step, which is 3.0e-12 1/s at this dt.
                        assert!(
                            carried_rate <= 16.0 * f64::EPSILON * density / dt,
                            "{model:?}, {gradient:?}: carried rate {carried_rate}"
                        );
                    }
                    assert!((stats.tau_range.0 - tau_at(carried_rate)).abs() < 4e-12);
                    assert!((stats.tau_range.1 - tau_at(carried_rate)).abs() < 4e-12);
                    assert_eq!((stats.floored, stats.capped), (0, 0));
                }
                // The two grids carry their own roundoff, so their taus agree to what the law makes of the
                // difference of the rates they carry (for a rotation), and to 4e-12 otherwise.
                assert!(
                    (planar.tau[0] - spatial.tau[1]).abs()
                        < 4e-12 + (tau_at(carried[0]) - tau_at(carried[1])).abs()
                );
            }
        }
    }

    #[test]
    fn dynamic_and_kinematic_si_inputs_apply_density_exactly_once() {
        let density = 800.0;
        let kinematic = [
            PhysicalRheology::Newtonian {
                viscosity: 0.08 / density,
            },
            PhysicalRheology::PowerLaw {
                consistency: 0.04 / density,
                n: 0.6,
                minimum_shear_rate_per_s: 0.002,
            },
            PhysicalRheology::Carreau {
                viscosity_zero: 0.4 / density,
                viscosity_infinite: 0.02 / density,
                relaxation_time_s: 0.03,
                n: 0.4,
            },
            PhysicalRheology::HerschelBulkley {
                consistency: 0.05 / density,
                n: 0.7,
                yield_stress: 0.2 / density,
                regularization_time_s: 0.6,
                minimum_shear_rate_per_s: 0.003,
            },
            PhysicalRheology::CarreauYasuda {
                viscosity_zero: 0.4 / density,
                viscosity_infinite: 0.02 / density,
                relaxation_time_s: 0.03,
                n: 0.4,
                a: 0.7,
            },
        ];
        for (dynamic, kinematic) in models().into_iter().zip(kinematic) {
            let dynamic = dynamic
                .convert(0.002, 0.002, density, LiquidViscosityBasis::Dynamic)
                .unwrap();
            let kinematic_same = kinematic
                .convert(0.002, 0.002, density, LiquidViscosityBasis::Kinematic)
                .unwrap();
            let kinematic_other = kinematic
                .convert(0.002, 0.002, 2.0 * density, LiquidViscosityBasis::Kinematic)
                .unwrap();
            for rate in [0.0, 1e-8, 0.002, 0.6] {
                let expected = dynamic.viscosity(rate).to_bits();
                assert_eq!(kinematic_same.viscosity(rate).to_bits(), expected);
                assert_eq!(kinematic_other.viscosity(rate).to_bits(), expected);
            }
        }
    }

    #[test]
    fn converted_coefficients_survive_cancelling_scales_but_refuse_lost_terms() {
        let nu = PhysicalRheology::Newtonian { viscosity: 1e-250 }
            .convert(1e-100, 1e-100, 1e100, LiquidViscosityBasis::Dynamic)
            .unwrap()
            .viscosity(0.0);
        assert!((nu / 1e-250 - 1.0).abs() < 1e-12);
        let power = PhysicalRheology::PowerLaw {
            consistency: 1e-250,
            n: 3.0,
            minimum_shear_rate_per_s: 1e-8,
        }
        .convert(1.0, 1e-100, 1.0, LiquidViscosityBasis::Kinematic)
        .unwrap();
        match power {
            Rheology::HerschelBulkley { k, .. } => assert!((k / 1e-150 - 1.0).abs() < 1e-12),
            _ => panic!("physical power law must preserve its explicit cutoff"),
        }
        let large = PhysicalRheology::Newtonian { viscosity: 1e200 }
            .convert(1e200, 1e200, 1.0, LiquidViscosityBasis::Kinematic)
            .unwrap();
        assert!((large.viscosity(0.0) - 1.0).abs() < 1e-14);
        for (viscosity, dx, dt) in [(1e308, 1e-200, 1e-200), (1e-300, 1e100, 1e-100)] {
            assert!(matches!(
                PhysicalRheology::Newtonian { viscosity }.convert(
                    dx,
                    dt,
                    1.0,
                    LiquidViscosityBasis::Kinematic
                ),
                Err(PhysicalRheologyError::Unrepresentable(_)),
            ));
        }
        assert!(matches!(
            PhysicalRheology::Carreau {
                viscosity_zero: 0.1,
                viscosity_infinite: 0.01,
                relaxation_time_s: 1e-300,
                n: 0.5,
            }
            .convert(1e50, 1e100, 1.0, LiquidViscosityBasis::Kinematic),
            Err(PhysicalRheologyError::Unrepresentable("constitutive time")),
        ));
        assert!(matches!(
            PhysicalRheology::PowerLaw {
                consistency: 0.1,
                n: 1.0,
                minimum_shear_rate_per_s: 1e-300,
            }
            .convert(1e-50, 1e-100, 1.0, LiquidViscosityBasis::Kinematic),
            Err(PhysicalRheologyError::Unrepresentable("minimum shear rate")),
        ));
        for invalid in [
            PhysicalRheology::Newtonian {
                viscosity: f64::INFINITY,
            },
            PhysicalRheology::PowerLaw {
                consistency: 0.1,
                n: 0.0,
                minimum_shear_rate_per_s: 1e-6,
            },
            PhysicalRheology::Carreau {
                viscosity_zero: 0.1,
                viscosity_infinite: 0.2,
                relaxation_time_s: 1.0,
                n: 2.0,
            },
            PhysicalRheology::HerschelBulkley {
                consistency: 0.1,
                n: 1.0,
                yield_stress: -1.0,
                regularization_time_s: 1.0,
                minimum_shear_rate_per_s: 1e-6,
            },
        ] {
            assert!(matches!(
                invalid.convert(0.01, 0.001, 800.0, LiquidViscosityBasis::Dynamic),
                Err(PhysicalRheologyError::Invalid(_))
            ));
        }
    }

    #[test]
    fn physical_binding_retains_source_scale_and_feeds_both_field_dimensions() {
        use fs_material::fluid_transport::FluidTransportState;
        let material = FluidTransportState::try_new(800.0, 2000.0, 0.2, 0.08).unwrap();
        let scale =
            crate::material::plan_scaling_from_material(&material, 0.002, 0.002, 0.01).unwrap();
        let bound = scale
            .bind_rheology(models()[0], LiquidViscosityBasis::Dynamic)
            .unwrap();
        assert!(core::ptr::eq(bound.scaling(), &scale));
        assert!(core::ptr::eq(bound.scaling().material(), &material));
        assert_eq!(bound.scaling().transport().density_kg_m3(), 800.0);
        assert_eq!(bound.basis(), LiquidViscosityBasis::Dynamic);
        assert_eq!(bound.lattice_law().viscosity(0.1), scale.flow().viscosity);
        assert_eq!(viscosity_si(bound.physical(), 10.0), 0.08);
        let laws = models().map(|model| {
            scale
                .bind_rheology(model, LiquidViscosityBasis::Dynamic)
                .unwrap()
                .lattice_law()
        });
        let field2 = super::super::RheologyField::new([1, 5], &laws).unwrap();
        let field3 = super::super::RheologyField::new([1, 1, 5], &laws).unwrap();
        for (i, law) in laws.iter().enumerate() {
            assert_eq!(field2.law(i).unwrap().viscosity(0.01), law.viscosity(0.01));
            assert_eq!(field3.law(i).unwrap().viscosity(0.01), law.viscosity(0.01));
        }
    }
}
