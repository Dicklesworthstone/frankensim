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

    fn models() -> [PhysicalRheology; 4] {
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
        let field2 = super::super::RheologyField::new([2, 2], &laws).unwrap();
        let field3 = super::super::RheologyField::new([1, 1, 4], &laws).unwrap();
        for (i, law) in laws.iter().enumerate() {
            assert_eq!(field2.law(i).unwrap().viscosity(0.01), law.viscosity(0.01));
            assert_eq!(field3.law(i).unwrap().viscosity(0.01), law.viscosity(0.01));
        }
    }
}
