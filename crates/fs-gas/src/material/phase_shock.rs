//! Homogeneous-equilibrium normal shocks with original phase enthalpy.
//!
//! Mass and momentum define the Rayleigh volume at each downstream pressure:
//! v2 = v1 - (p2-p1)/j^2. The original PH owner imposes the Hugoniot increment
//! Delta h = (p2-p1)*(v1+v2)/2. Matching its actual volume closes the jump;
//! the accepted state independently closes mass, momentum and stagnation h.
//! Species formation/entropy datums cancel in native finite increments.
//!
//! Conservation equations: NASA TN D-8274, "Flow Through Normal Shocks",
//! https://ntrs.nasa.gov/citations/19770008402. The caloric and phase closure
//! is the caller's original NASA/Peng--Robinson source, not that report's gas.
//!
//! This is a steady planar jump in the shock frame, with fixed total species
//! inventory, one phase velocity and complete thermal/phase equilibrium on
//! each side. Phase amounts may change or a phase may disappear across the
//! jump. The caller supplies a bracket on an ordinary compressive Hugoniot
//! branch; finite phase screening and endpoint sonic checks do not establish
//! global uniqueness or a shock-structure/relaxation model.

use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError,
    enthalpy::MAX_ENTHALPY_FLASH_EVALUATIONS,
    equilibrium::{
        PhaseEnthalpyControls, PhaseEnthalpySolution, ResolvedThermochemicalPhaseState,
        acoustics::{
            PHASE_ACOUSTIC_PT_EVALUATIONS, PhaseEquilibriumAcousticControls,
            PhaseEquilibriumAcoustics,
        },
        prepared::PreparedPhaseEnthalpyPath,
    },
};
use fs_qty::{Pressure, Temperature, Velocity};

pub const MAX_SHOCK_PRESSURE_EVALUATIONS: usize = 128;
/// One original inlet, two acoustic stencils, and bounded native PH queries.
pub const MAX_SHOCK_PT_EVALUATIONS: usize = 1
    + 2 * PHASE_ACOUSTIC_PT_EVALUATIONS
    + MAX_SHOCK_PRESSURE_EVALUATIONS * MAX_ENTHALPY_FLASH_EVALUATIONS;

/// Source domains, actual conservation tolerances and whole-call work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureShockControl {
    pub phase: PhaseEnthalpyControls,
    /// Must use exactly the PH phase policy. Both endpoint stencils must pass.
    pub acoustics: PhaseEquilibriumAcousticControls,
    /// Increasing pressures strictly above the inlet, excluding the identity.
    /// Both Rayleigh volumes must be positive; EOS-volume residual signs must
    /// run from nonpositive to nonnegative on the admitted compression branch.
    pub pressure_bracket: [Pressure; 2],
    pub absolute_volume_tolerance_m3_kg: f64,
    pub absolute_mass_flux_tolerance_kg_m2_s: f64,
    pub absolute_momentum_tolerance_pa: f64,
    pub absolute_energy_tolerance_j_kg: f64,
    /// Numerical allowance only; actual signed entropy change is retained.
    pub entropy_tolerance_j_kg_k: f64,
    /// Added to each absolute allowance using its physical transfer scale:
    /// inlet v, mass flux, pressure rise or enthalpy rise. Never absolute h.
    pub relative_tolerance: f64,
    pub maximum_pressure_evaluations: usize,
    pub maximum_pt_evaluations: usize,
}
impl MixtureShockControl {
    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2], pressure_bracket: [Pressure; 2]) -> Self {
        let acoustics = PhaseEquilibriumAcousticControls::default();
        let mut phase = PhaseEnthalpyControls::new(temperature_bracket);
        phase.enthalpy.flash_options = acoustics.phase.flash;
        phase.stability = acoustics.phase.stability;
        Self {
            phase,
            acoustics,
            pressure_bracket,
            absolute_volume_tolerance_m3_kg: 1e-12,
            absolute_mass_flux_tolerance_kg_m2_s: 1e-8,
            absolute_momentum_tolerance_pa: 1e-5,
            absolute_energy_tolerance_j_kg: 1e-5,
            entropy_tolerance_j_kg_k: 1e-7,
            relative_tolerance: 1e-10,
            maximum_pressure_evaluations: 64,
            maximum_pt_evaluations: MAX_SHOCK_PT_EVALUATIONS,
        }
    }

    pub fn validate(self) -> Result<(), MixtureShockError> {
        self.phase
            .validate()
            .map_err(|e| source_error(MixtureShockStage::Hugoniot, 0, e))?;
        self.acoustics
            .validate()
            .map_err(|e| source_error(MixtureShockStage::UpstreamAcoustics, 0, e))?;
        let [lo, hi] = self.pressure_bracket.map(|p| p.value());
        if self.phase.phase_options() != self.acoustics.phase
            || !positive(lo)
            || !positive(hi)
            || lo >= hi
            || [
                self.absolute_volume_tolerance_m3_kg,
                self.absolute_mass_flux_tolerance_kg_m2_s,
                self.absolute_momentum_tolerance_pa,
                self.absolute_energy_tolerance_j_kg,
                self.entropy_tolerance_j_kg_k,
            ]
            .into_iter()
            .any(|x| !positive(x))
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || !(2..=MAX_SHOCK_PRESSURE_EVALUATIONS).contains(&self.maximum_pressure_evaluations)
            || !(1..=MAX_SHOCK_PT_EVALUATIONS).contains(&self.maximum_pt_evaluations)
        {
            return Err(Error::Invalid(
                "matching phase policies, increasing pressures, positive tolerances and bounded shock work required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixtureShockStage {
    Inlet,
    Hugoniot,
    UpstreamAcoustics,
    DownstreamAcoustics,
}

/// No partially admitted jump is published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureShockError {
    Invalid(&'static str),
    Equilibrium {
        stage: MixtureShockStage,
        completed_pt_evaluations: usize,
        source: ThermochemicalFlashError,
    },
    WorkLimit {
        pressure_evaluations: usize,
        pt_evaluations: usize,
    },
    OutsideBracket {
        volume_residuals_m3_kg: [f64; 2],
    },
    NotConverged {
        pressure_bracket_pa: [f64; 2],
        volume_residual_m3_kg: f64,
    },
    Admissibility {
        upstream_mach: f64,
        downstream_mach: Option<f64>,
        entropy_change_j_kg_k: Option<f64>,
    },
    Conservation {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MixtureShockError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium mixture normal shock: {self:?}")
    }
}
impl std::error::Error for MixtureShockError {}
type Error = MixtureShockError;

/// Original phase endpoints and the actual physical jump defects.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureNormalShock {
    pub upstream: ResolvedThermochemicalPhaseState,
    /// Native original-inlet PH increment; its counters cover this query only.
    pub downstream: PhaseEnthalpySolution,
    pub upstream_acoustics: PhaseEquilibriumAcoustics,
    pub downstream_acoustics: PhaseEquilibriumAcoustics,
    pub upstream_speed: Velocity,
    pub downstream_speed: Velocity,
    pub mass_flux_kg_m2_s: f64,
    pub upstream_mach: f64,
    pub downstream_mach: f64,
    pub enthalpy_increase_j_kg: f64,
    pub entropy_change_j_kg_k: f64,
    pub specific_volume_residual_m3_kg: f64,
    pub mass_flux_residual_kg_m2_s: f64,
    pub momentum_flux_residual_pa: f64,
    pub hugoniot_residual_j_kg: f64,
    pub stagnation_enthalpy_residual_j_kg: f64,
    /// Last signed Rayleigh/EOS bracket, not an interval-certified bound.
    pub pressure_bracket_pa: [f64; 2],
    pub pressure_evaluations: usize,
    /// Whole call including one inlet preparation and both acoustic stencils.
    pub pt_evaluations: usize,
    pub eos_evaluation_upper_bound: usize,
}

/// Solve a classical compressive normal shock in the homogeneous-equilibrium limit.
///
/// `upstream_speed` is a positive magnitude relative to the shock. This query
/// neither advances a field nor chooses a laboratory shock velocity. All
/// kinetic-energy loss becomes the original total enthalpy rise, which already
/// contains pressure-flow work and latent energy; neither is added again.
///
/// The inlet must be supersonic and the accepted outlet subsonic according to
/// the native equilibrating sound speed. Entropy production, compression and
/// both independent energy equations are checked without clipping. Across the
/// discontinuity the phase count may change; differentiable endpoint acoustic
/// neighborhoods are still required. Nonclassical/undercompressive shocks,
/// nucleation, phase slip, finite-rate chemistry and shock thickness need
/// separate physical models. No ideal-gas or homogeneous fallback is used.
///
/// One inlet is prepared for the entire pressure root. PH queries reserve work
/// for the final acoustic stencil; only the two accepted endpoints need acoustic
/// derivatives. Refusal or cancellation returns no partial physical result.
pub fn equilibrium_mixture_normal_shock(
    model: &ThermochemicalPengRobinsonModel,
    upstream_temperature: Temperature,
    upstream_pressure: Pressure,
    upstream_speed: Velocity,
    mut controls: MixtureShockControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureNormalShock, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let p1 = upstream_pressure.value();
    let t1 = upstream_temperature.value();
    let u1 = upstream_speed.value();
    let [tlo, thi] = controls
        .phase
        .enthalpy
        .temperature_bracket
        .map(|t| t.value());
    let [plo, phi] = controls.pressure_bracket.map(|p| p.value());
    if !positive(p1) || !positive(t1) || !positive(u1) || t1 < tlo || t1 > thi || plo <= p1 {
        return Err(Error::Invalid(
            "positive inlet state/speed, admitted inlet T and strictly compressive pressure bracket required",
        ));
    }
    if controls.maximum_pt_evaluations < 1 + 2 * PHASE_ACOUSTIC_PT_EVALUATIONS {
        return Err(Error::WorkLimit {
            pressure_evaluations: 0,
            pt_evaluations: 0,
        });
    }
    controls.phase.enthalpy.enthalpy_tolerance_j_kg = controls
        .phase
        .enthalpy
        .enthalpy_tolerance_j_kg
        .min(0.125 * controls.absolute_energy_tolerance_j_kg);
    controls
        .phase
        .validate()
        .map_err(|e| source_error(MixtureShockStage::Hugoniot, 0, e))?;
    let path = model
        .prepare_phase_enthalpy_path(
            upstream_temperature,
            upstream_pressure,
            upstream_pressure,
            controls.phase,
            &mut cancelled,
        )
        .map_err(|e| source_error(MixtureShockStage::Inlet, 0, e))?;
    let mut work = Work {
        pressure: 0,
        pt: path.preparation_pt_evaluations(),
        eos: path.preparation_eos_evaluation_upper_bound(),
    };
    let initial = path
        .solve_increment(0.0, 0, &mut cancelled)
        .map_err(|e| source_error(MixtureShockStage::Inlet, work.pt, e))?;
    let upstream = initial.into_state();
    let upstream_acoustics = acoustic(
        model,
        &upstream,
        MixtureShockStage::UpstreamAcoustics,
        controls,
        &mut work,
        &mut cancelled,
    )?;
    let mach1 = u1 / upstream_acoustics.state().sound_speed_m_s();
    if !mach1.is_finite() || mach1 <= 1.0 {
        return Err(Error::Admissibility {
            upstream_mach: mach1,
            downstream_mach: None,
            entropy_change_j_kg_k: None,
        });
    }
    let v1 = 1.0 / upstream.density_kg_m3();
    let flux = upstream.density_kg_m3() * u1;
    if !positive(v1) || !positive(flux) || !positive(u1 - (phi - p1) / flux) {
        return Err(Error::Invalid(
            "shock mass flux and complete Rayleigh bracket must have positive finite volume",
        ));
    }
    let context = Context {
        path: &path,
        controls,
        p1,
        v1,
        u1,
        flux,
    };
    let mut lower = context.sample(plo, &mut work, &mut cancelled)?;
    let mut upper = context.sample(phi, &mut work, &mut cancelled)?;
    if lower.volume_residual > 0.0 || upper.volume_residual < 0.0 {
        return Err(Error::OutsideBracket {
            volume_residuals_m3_kg: [lower.volume_residual, upper.volume_residual],
        });
    }
    let (accepted, bracket) = loop {
        poll(&mut cancelled)?;
        let bracket = [lower.pressure(), upper.pressure()];
        let best = if lower.volume_residual.abs() <= upper.volume_residual.abs() {
            &lower
        } else {
            &upper
        };
        if context.balanced(best)? {
            break (
                if lower.volume_residual.abs() <= upper.volume_residual.abs() {
                    lower
                } else {
                    upper
                },
                bracket,
            );
        }
        // A central safeguarded secant uses native volume residuals, never a
        // fabricated heat capacity or phase derivative. Fall back to bisection
        // if arithmetic or an endpoint-hugging secant cannot contract the bracket.
        let width = bracket[1] - bracket[0];
        let fraction = -lower.volume_residual / (upper.volume_residual - lower.volume_residual);
        let p = bracket[0]
            + width
                * if (0.1..=0.9).contains(&fraction) {
                    fraction
                } else {
                    0.5
                };
        if !(p > bracket[0] && p < bracket[1]) {
            return Err(Error::NotConverged {
                pressure_bracket_pa: bracket,
                volume_residual_m3_kg: best.volume_residual,
            });
        }
        let sample = context.sample(p, &mut work, &mut cancelled)?;
        if sample.volume_residual <= 0.0 {
            lower = sample;
        } else {
            upper = sample;
        }
    };
    let downstream_acoustics = acoustic(
        model,
        accepted.outlet.state(),
        MixtureShockStage::DownstreamAcoustics,
        controls,
        &mut work,
        &mut cancelled,
    )?;
    let mach2 = accepted.speed / downstream_acoustics.state().sound_speed_m_s();
    if !positive(mach2)
        || mach2 >= 1.0
        || accepted.outlet.state().density_kg_m3() <= upstream.density_kg_m3()
        || accepted.enthalpy <= 0.0
        || accepted.entropy < -controls.entropy_tolerance_j_kg_k
    {
        return Err(Error::Admissibility {
            upstream_mach: mach1,
            downstream_mach: Some(mach2),
            entropy_change_j_kg_k: Some(accepted.entropy),
        });
    }
    context.enforce_balance(&accepted)?;
    poll(&mut cancelled)?;
    Ok(MixtureNormalShock {
        upstream,
        downstream: accepted.outlet,
        upstream_acoustics,
        downstream_acoustics,
        upstream_speed,
        downstream_speed: Velocity::new(accepted.speed),
        mass_flux_kg_m2_s: flux,
        upstream_mach: mach1,
        downstream_mach: mach2,
        enthalpy_increase_j_kg: accepted.enthalpy,
        entropy_change_j_kg_k: accepted.entropy,
        specific_volume_residual_m3_kg: accepted.volume_residual,
        mass_flux_residual_kg_m2_s: accepted.mass_residual,
        momentum_flux_residual_pa: accepted.momentum_residual,
        hugoniot_residual_j_kg: accepted.hugoniot_residual,
        stagnation_enthalpy_residual_j_kg: accepted.stagnation_residual,
        pressure_bracket_pa: bracket,
        pressure_evaluations: work.pressure,
        pt_evaluations: work.pt,
        eos_evaluation_upper_bound: work.eos,
    })
}

struct Work {
    pressure: usize,
    pt: usize,
    eos: usize,
}
impl Work {
    fn exhausted(&self) -> Error {
        Error::WorkLimit {
            pressure_evaluations: self.pressure,
            pt_evaluations: self.pt,
        }
    }
    fn charge(&mut self, pt: usize, eos: usize) -> Result<(), Error> {
        self.pt = self
            .pt
            .checked_add(pt)
            .ok_or(Error::Invalid("shock PT counter overflow"))?;
        self.eos = self
            .eos
            .checked_add(eos)
            .ok_or(Error::Invalid("shock EOS counter overflow"))?;
        Ok(())
    }
}
struct Sample {
    outlet: PhaseEnthalpySolution,
    speed: f64,
    enthalpy: f64,
    entropy: f64,
    volume_residual: f64,
    mass_residual: f64,
    momentum_residual: f64,
    hugoniot_residual: f64,
    stagnation_residual: f64,
}
impl Sample {
    fn pressure(&self) -> f64 {
        self.outlet.state().pressure_pa()
    }
}
struct Context<'m, 'p> {
    path: &'p PreparedPhaseEnthalpyPath<'m>,
    controls: MixtureShockControl,
    p1: f64,
    v1: f64,
    u1: f64,
    flux: f64,
}
impl Context<'_, '_> {
    fn sample(
        &self,
        pressure: f64,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sample, Error> {
        poll(cancelled)?;
        let remaining = self.controls.maximum_pt_evaluations.saturating_sub(work.pt);
        if work.pressure >= self.controls.maximum_pressure_evaluations
            || remaining <= PHASE_ACOUSTIC_PT_EVALUATIONS
        {
            return Err(work.exhausted());
        }
        let dp = pressure - self.p1;
        let drop = dp / self.flux;
        let rayleigh_speed = self.u1 - drop;
        let rayleigh_volume = rayleigh_speed / self.flux;
        // Equivalent to dp*(v1+vRayleigh)/2 without squaring the mass flux.
        let target = (0.5 * drop) * (self.u1 + rayleigh_speed);
        if !positive(dp)
            || !positive(drop)
            || !positive(rayleigh_speed)
            || !positive(rayleigh_volume)
            || !positive(target)
        {
            return Err(Error::Invalid(
                "unrepresentable positive Rayleigh/Hugoniot state",
            ));
        }
        work.pressure += 1;
        let outlet = self
            .path
            .solve_increment_at_pressure(
                Pressure::new(pressure),
                target,
                remaining - PHASE_ACOUSTIC_PT_EVALUATIONS,
                &mut *cancelled,
            )
            .map_err(|e| source_error(MixtureShockStage::Hugoniot, work.pt, e))?;
        work.charge(outlet.pt_evaluations(), outlet.eos_evaluation_upper_bound())?;
        let rho = outlet.state().density_kg_m3();
        let volume = 1.0 / rho;
        let speed = self.flux / rho;
        let actual_drop = self.u1 - speed;
        let enthalpy = outlet.accepted_increment_j_kg().ok_or(Error::Invalid(
            "shock requires original PH enthalpy increment",
        ))?;
        let entropy = outlet
            .entropy_change_j_kg_k()
            .ok_or(Error::Invalid("shock requires original entropy increment"))?;
        let sample = Sample {
            volume_residual: volume - rayleigh_volume,
            mass_residual: rho * speed - self.flux,
            momentum_residual: dp - self.flux * actual_drop,
            hugoniot_residual: enthalpy - (0.5 * dp) * (self.v1 + volume),
            stagnation_residual: enthalpy - (0.5 * actual_drop) * (self.u1 + speed),
            outlet,
            speed,
            enthalpy,
            entropy,
        };
        if !positive(volume)
            || !positive(speed)
            || [
                sample.volume_residual,
                sample.mass_residual,
                sample.momentum_residual,
                sample.hugoniot_residual,
                sample.stagnation_residual,
                enthalpy,
                entropy,
            ]
            .into_iter()
            .any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "nonfinite original shock state or conservation residual",
            ));
        }
        Ok(sample)
    }

    fn defects(&self, sample: &Sample) -> Result<[(&'static str, f64, f64); 5], Error> {
        let c = self.controls;
        let energy =
            c.absolute_energy_tolerance_j_kg + c.relative_tolerance * sample.enthalpy.abs();
        let result = [
            (
                "specific volume [m3/kg]",
                sample.volume_residual,
                c.absolute_volume_tolerance_m3_kg + c.relative_tolerance * self.v1,
            ),
            (
                "mass flux [kg/(m2 s)]",
                sample.mass_residual,
                c.absolute_mass_flux_tolerance_kg_m2_s + c.relative_tolerance * self.flux,
            ),
            (
                "momentum flux [Pa]",
                sample.momentum_residual,
                c.absolute_momentum_tolerance_pa
                    + c.relative_tolerance * (sample.pressure() - self.p1),
            ),
            ("Hugoniot enthalpy [J/kg]", sample.hugoniot_residual, energy),
            (
                "stagnation enthalpy [J/kg]",
                sample.stagnation_residual,
                energy,
            ),
        ];
        if result.iter().any(|(_, _, tolerance)| !positive(*tolerance)) {
            return Err(Error::Invalid("unrepresentable physical shock tolerance"));
        }
        Ok(result)
    }
    fn balanced(&self, sample: &Sample) -> Result<bool, Error> {
        Ok(self
            .defects(sample)?
            .into_iter()
            .all(|(_, r, tolerance)| r.abs() <= tolerance))
    }
    fn enforce_balance(&self, sample: &Sample) -> Result<(), Error> {
        for (quantity, residual, tolerance) in self.defects(sample)? {
            if residual.abs() > tolerance {
                return Err(Error::Conservation {
                    quantity,
                    residual,
                    tolerance,
                });
            }
        }
        Ok(())
    }
}

fn acoustic(
    model: &ThermochemicalPengRobinsonModel,
    state: &ResolvedThermochemicalPhaseState,
    stage: MixtureShockStage,
    controls: MixtureShockControl,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PhaseEquilibriumAcoustics, Error> {
    if controls.maximum_pt_evaluations.saturating_sub(work.pt) < PHASE_ACOUSTIC_PT_EVALUATIONS {
        return Err(work.exhausted());
    }
    let response = model
        .phase_equilibrium_acoustics(
            Temperature::new(state.temperature_k()),
            Pressure::new(state.pressure_pa()),
            controls.acoustics,
            &mut *cancelled,
        )
        .map_err(|e| source_error(stage, work.pt, e))?;
    work.charge(
        response.pt_evaluations(),
        response.eos_evaluation_upper_bound(),
    )?;
    if response.source() != state {
        return Err(Error::Invalid(
            "shock acoustics must retain the exact original caloric state",
        ));
    }
    Ok(response)
}
fn source_error(
    stage: MixtureShockStage,
    completed: usize,
    source: ThermochemicalFlashError,
) -> Error {
    match source {
        ThermochemicalFlashError::Cancelled => Error::Cancelled,
        source => Error::Equilibrium {
            stage,
            completed_pt_evaluations: completed,
            source,
        },
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
