//! Source-bound vibration, hereditary material loss and finite self-heating.
//!
//! A straight uniaxial specimen has the displacement field u(x)=x*q/L. Its
//! uniform strain is q/L, its internal force is A*sigma, and the Galerkin
//! inertia is rho*A*L/3 plus an explicitly supplied end payload. The SAME
//! source caloric density supplies both thermal mass and specimen inertia.
//! This is a one-coordinate axial reduction, not a resolved elastic wave model.
//!
//! Backward Euler condenses inertia to an affine endpoint stress constraint.
//! The existing Maxwell/enthalpy owner solves that constraint at EVERY thermal
//! trial, retaining the original branch history until the complete step passes.
//! Physical Maxwell dashpot work and signed reversible thermoelastic exchange
//! reach that same thermal owner. Endpoint numerical losses in material and
//! kinetic energy remain separate, nonthermal accounts.
//!
//! Source admission requires constant free strain/density, a single solid phase
//! and a globally affine spring-modulus scale. The original material owner
//! retains coefficient work and mechanical entropy when stiffness changes.
//! Relaxation times and heat capacity can also depend on temperature. No expansion,
//! phase change, spatial heat gradients, radiation of
//! sound or high-frequency structural accuracy is inferred from this reduction.

use core::{convert::Infallible, fmt};
use fs_exec::CancelGate;
use fs_material::phase::{
    ResolvedHeatCapacityCurve, UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::thermovisco::Energy;
use fs_material::thermovisco::source::{ResolvedThermoMaxwell, SourceError};
use fs_material::thermovisco::thermal::{
    CaloricPoint, FeedbackConfig, FeedbackError, MechanicalControl,
};

/// Fixed reference geometry and a distinct inertial end load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Geometry {
    pub length_m: f64,
    pub cross_section_m2: f64,
    /// Actual exposed area for the existing uniform thermal transport owner.
    pub thermal_surface_m2: f64,
    /// Lumped moving payload; its heat capacity is outside this specimen model.
    pub end_payload_kg: f64,
}

/// Independent numerical acceptance gates; temporal refinement is still needed.
#[derive(Clone, Copy, Debug)]
pub struct Controls {
    pub feedback: FeedbackConfig,
    pub force_tolerance_n: f64,
    /// Whole specimen+inertia energy residual allowance [J].
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum Error<E = Infallible> {
    Source(SourceError),
    Feedback(FeedbackError<E>),
    Invalid(&'static str),
    Cancelled,
    ForceResidual { residual_n: f64, tolerance_n: f64 },
    EnergyResidual { residual_j: f64, tolerance_j: f64 },
}
impl<E: fmt::Display> fmt::Display for Error<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Source(e) => write!(f, "thermovisco dynamics source: {e}"),
            Self::Feedback(e) => write!(f, "thermovisco dynamics: {e}"),
            Self::Invalid(s) => write!(f, "thermovisco dynamics: {s}"),
            Self::Cancelled => f.write_str("thermovisco dynamics: cancelled"),
            Self::ForceResidual {
                residual_n,
                tolerance_n,
            } => write!(
                f,
                "thermovisco dynamics force residual {residual_n} N exceeds {tolerance_n} N"
            ),
            Self::EnergyResidual {
                residual_j,
                tolerance_j,
            } => write!(
                f,
                "thermovisco dynamics energy residual {residual_j} J exceeds {tolerance_j} J"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for Error<E> {}

/// One accepted material/thermal/inertial trajectory. Sources and geometry
/// remain attached; replacing coefficients cannot reset or relabel its memory.
#[derive(Clone, Debug)]
pub struct Oscillator<'a> {
    mechanical_source: &'a ResolvedThermoMaxwell,
    caloric_source: &'a ResolvedHeatCapacityCurve,
    point: CaloricPoint<'a>,
    geometry: Geometry,
    effective_mass_kg: f64,
    velocity_m_s: f64,
}

/// Complete proposed endpoint, including the ORIGINAL thermal-owner report.
/// Assign `next` only after all coupled owners accept this same candidate.
#[derive(Debug)]
pub struct Step<'a, R> {
    pub next: Oscillator<'a>,
    /// Original Maxwell free-energy/work/loss audit [J/m3]. Its stored energy
    /// is Helmholtz free energy; `material_internal_energy_change_j` below
    /// includes the mechanical entropy contribution at changing temperature.
    pub material_energy: Energy,
    pub transport: UniformEnthalpyStep<R>,
    pub thermal_iterations: usize,
    /// Final nonnegative Maxwell dashpot loss [J], re-evaluated at accepted temperature.
    pub material_heat_j: f64,
    /// Signed reversible exchange into caloric storage [J], distinct from loss.
    pub reversible_material_heat_j: f64,
    /// Change of mechanical internal energy psi + T*s over the same step [J].
    pub material_internal_energy_change_j: f64,
    /// Change of the mechanical T*s contribution [J].
    pub mechanical_entropy_energy_change_j: f64,
    /// Total signed material heat supplied to the final transport trial [J].
    /// Feedback checks it against dashpot loss plus reversible exchange.
    pub applied_material_heat_j: f64,
    pub heat_defect_j: f64,
    pub thermal_residual_j: f64,
    pub kinetic_change_j: f64,
    pub inertial_numerical_loss_j: f64,
    /// Sum of material and inertial backward-Euler losses. Never material heat.
    pub numerical_loss_j: f64,
    /// Prescribed endpoint force times actual end displacement increment [J].
    pub external_work_j: f64,
    pub force_residual_n: f64,
    /// Delta kinetic + delta mechanical internal energy + delta enthalpy + numerical loss
    /// minus external endpoint work and actual boundary heat [J].
    pub combined_residual_j: f64,
}

impl<'a> Oscillator<'a> {
    /// Start a stress-free specimen with a declared velocity field x*v/L.
    /// `birth_strain` is its original stress-free reference, not a preloading
    /// operation. Subsequent forces excite vibration without resetting memory.
    /// Both models must retain the SAME original card and fixed source context.
    pub fn try_new(
        mechanical: &'a ResolvedThermoMaxwell,
        caloric: &'a ResolvedHeatCapacityCurve,
        geometry: Geometry,
        temperature_k: f64,
        birth_strain: f64,
        velocity_m_s: f64,
    ) -> Result<Self, Error> {
        if !positive(geometry.length_m)
            || !positive(geometry.cross_section_m2)
            || !positive(geometry.thermal_surface_m2)
            || !geometry.end_payload_kg.is_finite()
            || geometry.end_payload_kg < 0.0
            || !velocity_m_s.is_finite()
        {
            return Err(Error::Invalid(
                "positive finite geometry, nonnegative payload and finite velocity required",
            ));
        }
        let volume = finite(geometry.length_m * geometry.cross_section_m2)?;
        if volume <= 0.0 {
            return Err(Error::Invalid("specimen volume is not representable"));
        }
        let point = mechanical
            .caloric_point(
                caloric,
                temperature_k,
                birth_strain,
                volume,
                geometry.thermal_surface_m2,
            )
            .map_err(Error::Source)?;
        let effective_mass_kg = finite(point.mass_kg() / 3.0 + geometry.end_payload_kg)?;
        if effective_mass_kg <= 0.0
            || (velocity_m_s / geometry.length_m).abs() > point.max_strain_rate_s()
        {
            return Err(Error::Invalid(
                "nonpositive modal mass or initial velocity outside the material strain-rate envelope",
            ));
        }
        finite(0.5 * effective_mass_kg * velocity_m_s * velocity_m_s)?;
        Ok(Self {
            mechanical_source: mechanical,
            caloric_source: caloric,
            point,
            geometry,
            effective_mass_kg,
            velocity_m_s,
        })
    }

    #[must_use]
    pub const fn point(&self) -> &CaloricPoint<'a> {
        &self.point
    }
    #[must_use]
    pub const fn mechanical_source(&self) -> &'a ResolvedThermoMaxwell {
        self.mechanical_source
    }
    #[must_use]
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.caloric_source
    }
    #[must_use]
    pub const fn geometry(&self) -> Geometry {
        self.geometry
    }
    /// Integral rho*(x/L)^2 dV plus the declared end payload [kg].
    #[must_use]
    pub const fn effective_mass_kg(&self) -> f64 {
        self.effective_mass_kg
    }
    #[must_use]
    pub const fn velocity_m_s(&self) -> f64 {
        self.velocity_m_s
    }
    #[must_use]
    pub fn extension_m(&self) -> f64 {
        self.geometry.length_m * self.point.mechanical().total_strain()
    }
    #[must_use]
    pub fn kinetic_energy_j(&self) -> f64 {
        0.5 * self.effective_mass_kg * self.velocity_m_s * self.velocity_m_s
    }

    /// Propose one force-driven step with simultaneous thermal/material feedback.
    ///
    /// m*(v1-v0)/dt + A*sigma1 = F1 and L*(epsilon1-epsilon0)=dt*v1
    /// give the affine load line used by the ORIGINAL Maxwell solver:
    /// sigma1 = F1/A + m*v0/(A*dt) - m*L/(A*dt^2)*(epsilon1-epsilon0).
    /// The old strain, branch history, temperature and duration are unchanged
    /// through all trials. Re-evaluating that line at the accepted temperature
    /// closes inertia and relaxation without stale-temperature mechanical forces.
    ///
    /// Endpoint work obeys F1*delta_q = delta_K + m*(delta_v)^2/2
    /// + material_work. The nonnegative kinetic term and the original Maxwell
    /// numerical term never enter transport. The material thermal owner supplies
    /// physical dashpot heat plus signed reversible exchange from its affine
    /// free-energy law; its mechanical internal energy closes the balance.
    ///
    /// The transport callback follows the existing candidate contract: it must
    /// not mutate reservoirs during trials. For example use the existing
    /// `LumpedThermalEnvironment::advance(cx, input)` and its Biot admission.
    /// Failure/cancellation returns no candidate and never changes this state.
    pub fn step<R, E>(
        &self,
        external_force_n: f64,
        duration_s: f64,
        controls: Controls,
        transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        gate: &CancelGate,
    ) -> Result<Step<'a, R>, Error<E>> {
        poll(gate)?;
        if !external_force_n.is_finite()
            || !positive(duration_s)
            || !positive(controls.force_tolerance_n)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(Error::Invalid(
                "finite force and positive interval/force/energy gates required",
            ));
        }
        let area = self.geometry.cross_section_m2;
        let inertia_per_area_dt = finite(self.effective_mass_kg / area / duration_s)?;
        let stiffness_pa = finite(inertia_per_area_dt * self.geometry.length_m / duration_s)?;
        let stress_at_reference_pa =
            finite(external_force_n / area + finite(inertia_per_area_dt * self.velocity_m_s)?)?;
        let stress_tolerance = finite(controls.force_tolerance_n / area)?;
        if stiffness_pa <= 0.0 || stress_tolerance <= 0.0 {
            return Err(Error::Invalid(
                "condensed inertia or force tolerance underflowed",
            ));
        }
        let old_strain = self.point.mechanical().total_strain();
        let feedback = self
            .point
            .advance_controlled(
                MechanicalControl::AffineStress {
                    reference_strain: old_strain,
                    stress_at_reference_pa,
                    stiffness_pa,
                    tolerance_pa: stress_tolerance,
                },
                duration_s,
                controls.feedback,
                transport,
                || gate.is_requested(),
            )
            .map_err(Error::Feedback)?;
        let displacement = finite(
            self.geometry.length_m
                * finite(feedback.next.mechanical().total_strain() - old_strain)?,
        )?;
        let velocity = finite(displacement / duration_s)?;
        finite(0.5 * self.effective_mass_kg * velocity * velocity)?;
        let delta_velocity = finite(velocity - self.velocity_m_s)?;
        let force_residual_n = finite(
            finite(self.effective_mass_kg * delta_velocity / duration_s)?
                + finite(area * feedback.next.mechanical().stress_pa())?
                - external_force_n,
        )?;
        if force_residual_n.abs() > controls.force_tolerance_n {
            return Err(Error::ForceResidual {
                residual_n: force_residual_n,
                tolerance_n: controls.force_tolerance_n,
            });
        }
        let volume = self.point.volume_m3();
        let kinetic_change_j = finite(
            0.5 * self.effective_mass_kg * delta_velocity * finite(velocity + self.velocity_m_s)?,
        )?;
        let inertial_numerical_loss_j =
            finite(0.5 * self.effective_mass_kg * delta_velocity * delta_velocity)?;
        let numerical_loss_j = finite(
            inertial_numerical_loss_j
                + finite(volume * feedback.mechanical_energy.numerical_dissipation)?,
        )?;
        let external_work_j = finite(external_force_n * displacement)?;
        let thermal_change = finite(
            self.point.mass_kg()
                * finite(
                    feedback.next.caloric_state().specific_enthalpy_j_kg()
                        - self.point.caloric_state().specific_enthalpy_j_kg(),
                )?,
        )?;
        let material_heat_j = finite(volume * feedback.mechanical_energy.dissipated)?;
        let combined_residual_j = finite(
            kinetic_change_j
                + feedback.mechanical_internal_energy_change_j
                + thermal_change
                + numerical_loss_j
                - external_work_j
                - feedback.transport.external_heat_j,
        )?;
        if combined_residual_j.abs() > controls.energy_tolerance_j {
            return Err(Error::EnergyResidual {
                residual_j: combined_residual_j,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        poll(gate)?;
        let next = Self {
            mechanical_source: self.mechanical_source,
            caloric_source: self.caloric_source,
            point: feedback.next,
            geometry: self.geometry,
            effective_mass_kg: self.effective_mass_kg,
            velocity_m_s: velocity,
        };
        Ok(Step {
            next,
            material_energy: feedback.mechanical_energy,
            transport: feedback.transport,
            thermal_iterations: feedback.iterations,
            material_heat_j,
            reversible_material_heat_j: feedback.reversible_heat_j,
            material_internal_energy_change_j: feedback.mechanical_internal_energy_change_j,
            mechanical_entropy_energy_change_j: feedback.mechanical_entropy_energy_change_j,
            applied_material_heat_j: feedback.applied_heat_j,
            heat_defect_j: feedback.heat_defect_j,
            thermal_residual_j: feedback.thermal_residual_j,
            kinetic_change_j,
            inertial_numerical_loss_j,
            numerical_loss_j,
            external_work_j,
            force_residual_n,
            combined_residual_j,
        })
    }
}

fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn finite<E>(v: f64) -> Result<f64, Error<E>> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(Error::Invalid(
            "unrepresentable inertial or energy arithmetic",
        ))
    }
}
fn poll<E>(gate: &CancelGate) -> Result<(), Error<E>> {
    if gate.is_requested() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
