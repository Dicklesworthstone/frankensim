//! Rate-independent dielectric switching through the shared scalar play solver.
//!
//! `D = epsilon_b E + sum(P_k)`,
//! `W = epsilon_b E²/2 + sum(P_k²/(2 chi_k))`, and
//! `Q = sum(Ec_k |delta P_k|)`. Each polarization minimizes
//! `P²/(2 chi) - E P + Ec |P-P_old|` on `[-Ps, Ps]`.
//! The numerical problem is exactly the existing scalar magnetic play problem;
//! only its SI quantities and terminal reduction differ. Its return map,
//! safeguarded inverse, exact segment work and history admission are reused.
//! Magnetic types never escape this electric interface.
//!
//! This quadratic-potential, bounded-polarization model follows the dry-friction
//! energy construction of Babori, Barati & Daniel (2024),
//! doi:10.1016/j.euromechsol.2023.105151, with an explicitly simpler reversible
//! potential. It does not reproduce their multiscale anhysteretic law or confer
//! measured ferroelectric accuracy on supplied fits. Temperature, geometry and
//! stress stay fixed; dielectric breakdown, fatigue, leakage, rate dependence,
//! mechanical switching strain and Curie transitions require separate owners.
//!
//! Supply every field/displacement reversal as a separate monotone segment.
//! The exact integral `E dD = delta W + Q` includes pinning and saturation
//! crossings. A trapezoid between endpoint fields generally does not. This
//! history-dependent law is not a path-independent storage function of D.

use crate::field_storage::{
    FieldStorageError,
    hysteresis::{HysteresisBranch, HysteresisStep, MagneticHysteresisState, ScalarHysteresis},
};

/// Original material-card roles and source-bound fitted coefficients.
pub mod source;

/// One explicitly fitted polarization contribution, with no hidden weights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricHysteresisBranch {
    /// Reversible dP/dE while this branch switches [F/m], not relative epsilon.
    pub susceptibility_f_m: f64,
    /// Nonnegative branch pinning field [V/m], not an entire loop's coercivity.
    pub pinning_field_v_m: f64,
    /// Positive bound on this branch's |P| [C/m²].
    pub saturation_polarization_c_m2: f64,
}

/// Immutable fixed-condition uniaxial law with positive background permittivity.
#[derive(Clone, Debug, PartialEq)]
pub struct ScalarDielectricHysteresis {
    play: ScalarHysteresis,
}

/// Accepted electric field, displacement and complete switching memory.
/// Private fields prevent history transplantation between different fits.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricHysteresisState {
    play: MagneticHysteresisState,
}

impl DielectricHysteresisState {
    pub const fn field_v_m(&self) -> f64 {
        self.play.field_a_m()
    }
    pub const fn displacement_c_m2(&self) -> f64 {
        self.play.flux_density_t()
    }
    pub fn polarizations_c_m2(&self) -> &[f64] {
        self.play.polarizations_t()
    }
    pub const fn stored_energy_j_m3(&self) -> f64 {
        self.play.stored_energy_j_m3()
    }
}

/// Pure trial response. Retaining `state` commits the complete material history.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricHysteresisStep {
    pub state: DielectricHysteresisState,
    /// Physical switching loss, deposited at most once [J/m³].
    pub irreversible_heat_j_m3: f64,
    /// Exact integral E dD on this monotone segment [J/m³].
    pub electric_work_j_m3: f64,
    /// Stable change of original reversible storage [J/m³].
    pub energy_change_j_m3: f64,
    /// Work divided by actual delta D; endpoint E at zero increment [V/m].
    pub mean_field_v_m: f64,
    /// Frozen-history dD/dE [F/m]. Kinks select the pinned/saturated side.
    pub incremental_permittivity_f_m: f64,
}

impl DielectricHysteresisStep {
    fn from_play(step: HysteresisStep) -> Self {
        Self {
            state: DielectricHysteresisState { play: step.state },
            irreversible_heat_j_m3: step.irreversible_heat_j_m3,
            electric_work_j_m3: step.magnetic_work_j_m3,
            energy_change_j_m3: step.energy_change_j_m3,
            mean_field_v_m: step.mean_field_a_m,
            incremental_permittivity_f_m: step.incremental_permeability_h_m,
        }
    }
}

impl ScalarDielectricHysteresis {
    /// Construct from explicit SI coefficients. At most 1024 branches; inverse
    /// work is bounded by the shared 96-iteration safeguarded active-set solve.
    /// Zero branches is the ordinary linear, lossless dielectric limit.
    pub fn new(
        background_permittivity_f_m: f64,
        field_limit_v_m: f64,
        branches: Vec<DielectricHysteresisBranch>,
    ) -> Result<Self, FieldStorageError> {
        if branches.len() > 1024 {
            return Err(FieldStorageError::Invalid(
                "dielectric switching exceeds native branch capacity",
            ));
        }
        Ok(Self {
            play: ScalarHysteresis::new(
                background_permittivity_f_m,
                field_limit_v_m,
                branches
                    .into_iter()
                    .map(|b| HysteresisBranch {
                        susceptibility_h_m: b.susceptibility_f_m,
                        pinning_field_a_m: b.pinning_field_v_m,
                        saturation_polarization_t: b.saturation_polarization_c_m2,
                    })
                    .collect(),
            )?,
        })
    }

    /// Virgin unpolarized state. Prepare poled specimens through their history.
    pub fn initial_state(&self) -> DielectricHysteresisState {
        DielectricHysteresisState {
            play: self.play.initial_state(),
        }
    }

    /// Exact field-driven monotone segment, including remanence and minor loops.
    pub fn candidate_field(
        &self,
        previous: &DielectricHysteresisState,
        field_v_m: f64,
    ) -> Result<DielectricHysteresisStep, FieldStorageError> {
        self.play
            .candidate_field(&previous.play, field_v_m)
            .map(DielectricHysteresisStep::from_play)
    }

    /// Invert the SAME frozen-history law; only the final branch state allocates.
    /// Returned D is the actual constitutive response, never a clipped target.
    pub fn candidate_displacement(
        &self,
        previous: &DielectricHysteresisState,
        displacement_c_m2: f64,
    ) -> Result<DielectricHysteresisStep, FieldStorageError> {
        self.play
            .candidate_flux_density(&previous.play, displacement_c_m2)
            .map(DielectricHysteresisStep::from_play)
    }
}

/// Uniform fixed parallel-plate capacitor: q=A D, V=d E, and volume=A d.
/// Electrode/free charge includes the background and every polarization branch.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformHystereticCapacitor {
    law: ScalarDielectricHysteresis,
    source: Option<source::ResolvedDielectricHysteresis>,
    area_m2: f64,
    gap_m: f64,
    volume_m3: f64,
}

/// Terminal response with physical switching heat distinct from stored energy.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticCapacitorStep {
    pub material: DielectricHysteresisStep,
    pub charge_c: f64,
    pub voltage_v: f64,
    /// Integral V dq divided by delta q, not the endpoint voltage [V].
    pub mean_voltage_v: f64,
    pub incremental_capacitance_f: f64,
    pub electric_work_j: f64,
    pub switching_heat_j: f64,
    pub stored_energy_j: f64,
}

impl UniformHystereticCapacitor {
    pub fn new(
        law: ScalarDielectricHysteresis,
        area_m2: f64,
        gap_m: f64,
    ) -> Result<Self, FieldStorageError> {
        if [area_m2, gap_m, area_m2 * gap_m, area_m2 / gap_m]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(FieldStorageError::Invalid(
                "hysteretic capacitor requires positive representable area, gap and volume",
            ));
        }
        let capacitor = Self {
            law,
            source: None,
            area_m2,
            gap_m,
            volume_m3: area_m2 * gap_m,
        };
        capacitor.hold_candidate(&capacitor.initial_state())?;
        Ok(capacitor)
    }

    pub fn initial_state(&self) -> DielectricHysteresisState {
        self.law.initial_state()
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    pub const fn law(&self) -> &ScalarDielectricHysteresis {
        &self.law
    }

    /// Admit the same memory and recover terminals without switching or heat.
    pub fn hold_candidate(
        &self,
        state: &DielectricHysteresisState,
    ) -> Result<HystereticCapacitorStep, FieldStorageError> {
        self.finish(self.law.candidate_field(state, state.field_v_m())?)
    }

    pub fn candidate_voltage(
        &self,
        previous: &DielectricHysteresisState,
        voltage_v: f64,
    ) -> Result<HystereticCapacitorStep, FieldStorageError> {
        self.finish(self.law.candidate_field(previous, voltage_v / self.gap_m)?)
    }

    /// Charge-controlled trial for a native circuit. Every nonlinear iteration
    /// must reuse the original accepted `previous`, never another trial's state.
    pub fn candidate_charge(
        &self,
        previous: &DielectricHysteresisState,
        charge_c: f64,
    ) -> Result<HystereticCapacitorStep, FieldStorageError> {
        self.finish(
            self.law
                .candidate_displacement(previous, charge_c / self.area_m2)?,
        )
    }

    fn finish(
        &self,
        material: DielectricHysteresisStep,
    ) -> Result<HystereticCapacitorStep, FieldStorageError> {
        let result = HystereticCapacitorStep {
            charge_c: self.area_m2 * material.state.displacement_c_m2(),
            voltage_v: self.gap_m * material.state.field_v_m(),
            mean_voltage_v: self.gap_m * material.mean_field_v_m,
            incremental_capacitance_f: self.area_m2 / self.gap_m
                * material.incremental_permittivity_f_m,
            electric_work_j: self.volume_m3 * material.electric_work_j_m3,
            switching_heat_j: self.volume_m3 * material.irreversible_heat_j_m3,
            stored_energy_j: self.volume_m3 * material.state.stored_energy_j_m3(),
            material,
        };
        if [
            result.charge_c,
            result.voltage_v,
            result.mean_voltage_v,
            result.incremental_capacitance_f,
            result.electric_work_j,
            result.switching_heat_j,
            result.stored_energy_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || result.incremental_capacitance_f <= 0.0
        {
            return Err(FieldStorageError::Invalid(
                "hysteretic capacitor response is unrepresentable",
            ));
        }
        Ok(result)
    }
}
