//! Finite-rate ion channels beside independent, fully selective water pores.
//!
//! Fluxes are positive left to right; `psi=F*(phiL-phiR)/(RT)`. Each channel
//! integrates `J_i=P_i*(A_i*c_i-dc_i/ds)` on the unit membrane coordinate,
//! with `A_i=ln(gammaL_i/gammaR_i)+z_i*psi` constant along that coordinate.
//! Thus `J_i=P_i*[B(-A_i)cL_i-B(A_i)cR_i]`, `B(x)=x/(exp(x)-1)`.
//! The effective electric AND excess chemical potentials are linear through
//! the membrane. This is not an exact solution of a nonlinear membrane interior.
//! Bulk ideal or virial activities and osmotic pressure retain their original
//! free-energy owner. Permeances include the caller's membrane thickness and
//! partition factors; no diffusivity, dissociation or extra standard state is
//! inferred. The constant-field precedent is Goldman (1943),
//! <https://doi.org/10.1085/jgp.27.1.37>.
//!
//! The floating voltage solves `sum(z_i*J_i)=0` with a monotone bounded root.
//! Ionic fluxes share their conjugate electrochemical forces, so each
//! `RT*J_i*[ln(cL_i/cR_i)+A_i]` is nonnegative. An exact one-sided-zero
//! concentration has finite flux but a singular instantaneous dissipation
//! diagnostic; no concentration floor supplies its logarithm. Exact blocked
//! counterion-only limits retain infinite potential and zero charged flux.
//!
//! Water uses `Jv=Lp*(Delta(p)-Delta(Pi))` through separate selective pores.
//! These stationary-solvent ion channels add no solvent drag or inferred
//! reciprocal coupling. The model is isothermal, neglects solute partial
//! volumes and stores no membrane charge, heat or pore inventory. Finite bulk
//! electroneutrality, fixed charge and conservative time integration belong to
//! the consumer. An isolated pair may use this local open-circuit closure;
//! arbitrary conducting cycles require globally compatible electric potentials.

use super::donnan::{DonnanPotential, FARADAY_C_PER_MOL};
use super::nonideal::log_ratio;
use super::{
    AqueousSolution, HydraulicPermeance, MAX_OSMOTIC_SOLUTES, OsmosisError, SolutePermeance, Sum,
    det, diffusive_exchange, finite, layout, nonnegative, positive, shape, thermal_scale,
};

#[cfg(test)]
mod tests;

/// Original thermodynamic admission, bounded current solve, or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum ElectrodiffusionError {
    Material(OsmosisError),
    Invalid(&'static str),
    NonConvergence,
    Cancelled,
}
impl core::fmt::Display for ElectrodiffusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electrodiffusive membrane: {self:?}")
    }
}
impl std::error::Error for ElectrodiffusionError {}
impl From<OsmosisError> for ElectrodiffusionError {
    fn from(value: OsmosisError) -> Self {
        Self::Material(value)
    }
}
type Result<T> = core::result::Result<T, ElectrodiffusionError>;

/// Accuracy of actual charge flux, independently of any inventory/time budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrodiffusionControls {
    /// Positive absolute maximum elementary-charge-equivalent flux [mol/(m² s)].
    pub absolute_current_tolerance_mol_m2_s: f64,
    /// Both this relative test and the absolute test must pass. Its scale is
    /// `sum |z_i| P_i max(cL_i,cR_i)`, independent of the trial voltage.
    pub relative_current_tolerance: f64,
    /// All candidate current/flux evaluations, including rejected root trials.
    pub maximum_potential_evaluations: usize,
}
impl Default for ElectrodiffusionControls {
    fn default() -> Self {
        Self {
            absolute_current_tolerance_mol_m2_s: 1e-14,
            relative_current_tolerance: 4e-13,
            maximum_potential_evaluations: 128,
        }
    }
}

/// Common original aqueous model and an already reference-cancelled load.
#[derive(Clone, Copy, Debug)]
pub struct ElectrodiffusionProblem<'a> {
    pub left: &'a AqueousSolution,
    pub right: &'a AqueousSolution,
    /// Hydrostatic `p_left-p_right` [Pa], with the common gauge cancelled first.
    pub pressure_drop_pa: f64,
    /// Additional ideal retained particles [mol], not already in either solution.
    /// Fixed charge is not itself an additional osmotic particle count.
    pub trapped_particles_mol: [f64; 2],
}

/// Actual open-circuit constitutive result; all fluxes are per membrane area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrodiffusionReport {
    pub volume_flux_m_s: f64,
    pub osmotic_pressure_drop_pa: f64,
    pub driving_pressure_pa: f64,
    /// LEFT-minus-RIGHT potential; an absent conducting pool is indeterminate.
    pub potential: DonnanPotential,
    pub charge_flux_residual_mol_m2_s: f64,
    pub potential_evaluations: usize,
    pub hydraulic_dissipation_w_m2: f64,
    /// Electrochemical dissipation; None for a conducting one-sided-zero front.
    pub ionic_dissipation_w_m2: Option<f64>,
    pub total_dissipation_w_m2: Option<f64>,
}

/// Constant explicit permeances on the original, at-most-64-solute ordering.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectrodiffusiveMembrane {
    hydraulic: HydraulicPermeance,
    valences: Box<[i32]>,
    permeances: Box<[SolutePermeance]>,
}
struct Face {
    left: [f64; MAX_OSMOTIC_SOLUTES],
    right: [f64; MAX_OSMOTIC_SOLUTES],
    excess_drop: [f64; MAX_OSMOTIC_SOLUTES],
    rt: f64,
    osmotic_drop: f64,
    current_scale: f64,
    positive_carrier: bool,
    negative_carrier: bool,
}
#[derive(Clone, Copy)]
struct CurrentPoint {
    psi: f64,
    current: f64,
    derivative: f64,
}

impl ElectrodiffusiveMembrane {
    /// Zero permeance excludes an ion channel. Zero valence permits an explicit
    /// neutral cosolute. Empty/all-blocked conducting layouts remain meaningful
    /// water membranes; they do not determine an electric potential.
    pub fn new(
        hydraulic: HydraulicPermeance,
        valences: Vec<i32>,
        permeances: Vec<SolutePermeance>,
    ) -> Result<Self> {
        layout(valences.len())?;
        shape(valences.len(), permeances.len())?;
        Ok(Self {
            hydraulic,
            valences: valences.into_boxed_slice(),
            permeances: permeances.into_boxed_slice(),
        })
    }
    pub fn valences(&self) -> &[i32] {
        &self.valences
    }
    pub fn permeances(&self) -> &[SolutePermeance] {
        &self.permeances
    }
    pub const fn hydraulic_permeance(&self) -> HydraulicPermeance {
        self.hydraulic
    }
    pub fn heap_bytes(&self) -> usize {
        core::mem::size_of_val(self.valences.as_ref())
            + core::mem::size_of_val(self.permeances.as_ref())
    }

    /// Allocation-free evaluation with atomic publication of the output slice.
    pub fn flux_into(
        &self,
        problem: ElectrodiffusionProblem<'_>,
        controls: ElectrodiffusionControls,
        output: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ElectrodiffusionReport> {
        self.flux_into_counted(problem, controls, output, &mut 0, cancelled)
    }

    /// The external counter also records trials from unsuccessful calls. The
    /// per-call cap includes every attempted current/flux evaluation, including
    /// a closed-form absent-pool or infinite-potential boundary evaluation.
    /// Scratch and counters may change on refusal; caller fluxes never do.
    pub fn flux_into_counted(
        &self,
        problem: ElectrodiffusionProblem<'_>,
        controls: ElectrodiffusionControls,
        output: &mut [f64],
        potential_evaluations: &mut usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ElectrodiffusionReport> {
        shape(self.valences.len(), output.len())?;
        if !positive(controls.absolute_current_tolerance_mol_m2_s)
            || !positive(controls.relative_current_tolerance)
            || controls.relative_current_tolerance >= 1.0
            || controls.maximum_potential_evaluations == 0
            || controls.maximum_potential_evaluations > 4096
        {
            return Err(ElectrodiffusionError::Invalid(
                "explicit positive current accuracy and bounded potential work required",
            ));
        }
        poll(&mut cancelled)?;
        let face = self.prepare(problem)?;
        let tolerance = controls
            .absolute_current_tolerance_mol_m2_s
            .min(controls.relative_current_tolerance * face.current_scale);
        let mut candidate = [0.0; MAX_OSMOTIC_SOLUTES];
        let mut evaluations = 0;
        let mut evaluate = |psi: f64| -> Result<CurrentPoint> {
            poll(&mut cancelled)?;
            if evaluations >= controls.maximum_potential_evaluations {
                return Err(ElectrodiffusionError::NonConvergence);
            }
            evaluations += 1;
            *potential_evaluations =
                potential_evaluations
                    .checked_add(1)
                    .ok_or(ElectrodiffusionError::Invalid(
                        "potential work counter overflow",
                    ))?;
            self.current_at(&face, psi, &mut candidate)
        };
        let (potential, current) = match (face.positive_carrier, face.negative_carrier) {
            (true, true) => {
                let scale = self
                    .valences
                    .iter()
                    .map(|z| f64::from(*z).abs())
                    .fold(1.0, f64::max);
                let point = solve_current(tolerance, 1.0 / scale, &mut evaluate)?;
                (
                    DonnanPotential::Finite {
                        dimensionless_left_minus_right: point.psi,
                        voltage_left_minus_right_v: finite(
                            face.rt * point.psi / FARADAY_C_PER_MOL,
                            "membrane voltage is not representable",
                        )?,
                    },
                    point.current,
                )
            }
            (positive, negative) => {
                let potential = if positive {
                    DonnanPotential::NegativeInfinite
                } else if negative {
                    DonnanPotential::PositiveInfinite
                } else {
                    DonnanPotential::Indeterminate
                };
                // An infinite field blocks all conducting charged species in
                // this one-direction-only limit. Neutral channels still flow.
                let psi = if positive {
                    f64::NEG_INFINITY
                } else if negative {
                    f64::INFINITY
                } else {
                    0.0
                };
                let point = evaluate(psi)?;
                (potential, point.current)
            }
        };
        let driving_pressure_pa = finite(
            problem.pressure_drop_pa - face.osmotic_drop,
            "membrane driving pressure is not representable",
        )?;
        let volume_flux_m_s = finite(
            self.hydraulic.value() * driving_pressure_pa,
            "membrane water flux is not representable",
        )?;
        let hydraulic_dissipation_w_m2 = finite(
            volume_flux_m_s * driving_pressure_pa,
            "hydraulic membrane dissipation is not representable",
        )?;
        let ionic_dissipation_w_m2 = self.dissipation(&face, potential, &candidate)?;
        let total_dissipation_w_m2 = ionic_dissipation_w_m2
            .map(|ionic| {
                finite(
                    hydraulic_dissipation_w_m2 + ionic,
                    "total membrane dissipation is not representable",
                )
            })
            .transpose()?;
        let report = ElectrodiffusionReport {
            volume_flux_m_s,
            osmotic_pressure_drop_pa: face.osmotic_drop,
            driving_pressure_pa,
            potential,
            charge_flux_residual_mol_m2_s: current,
            potential_evaluations: evaluations,
            hydraulic_dissipation_w_m2,
            ionic_dissipation_w_m2,
            total_dissipation_w_m2,
        };
        poll(&mut cancelled)?;
        output.copy_from_slice(&candidate[..output.len()]);
        Ok(report)
    }

    fn prepare(&self, problem: ElectrodiffusionProblem<'_>) -> Result<Face> {
        let n = self.valences.len();
        for solution in [problem.left, problem.right] {
            shape(n, solution.solute_moles().len())?;
        }
        if problem.left.temperature_k() != problem.right.temperature_k()
            || problem.left.solvent_molar_volume_m3_mol()
                != problem.right.solvent_molar_volume_m3_mol()
            || problem.left.nonideal_model() != problem.right.nonideal_model()
            || !problem.pressure_drop_pa.is_finite()
            || problem
                .trapped_particles_mol
                .iter()
                .any(|n| !nonnegative(*n))
        {
            return Err(ElectrodiffusionError::Invalid(
                "common original aqueous model, finite pressure drop and nonnegative trapped osmoles required",
            ));
        }
        let mut face = Face {
            left: [0.0; MAX_OSMOTIC_SOLUTES],
            right: [0.0; MAX_OSMOTIC_SOLUTES],
            excess_drop: [0.0; MAX_OSMOTIC_SOLUTES],
            rt: thermal_scale(problem.left.temperature_k())?,
            osmotic_drop: 0.0,
            current_scale: 0.0,
            positive_carrier: false,
            negative_carrier: false,
        };
        problem.left.concentrations_into(&mut face.left[..n])?;
        problem.right.concentrations_into(&mut face.right[..n])?;
        let mut pressure = Sum::default();
        let mut scale = Sum::default();
        for i in 0..n {
            pressure.add(face.left[i] - face.right[i]);
            if let Some(model) = problem.left.nonideal_model() {
                // B*(cL-cR) cancels common excess references BEFORE arithmetic.
                // It is the exact difference from this original free energy,
                // including for absent species whose ideal logarithm is singular.
                let mut excess = Sum::default();
                for j in 0..n {
                    excess.add(
                        model.interaction_matrix_m3_mol()[i * n + j]
                            * (face.left[j] - face.right[j]),
                    );
                }
                face.excess_drop[i] = excess.finish()?;
                // Symmetry gives Delta(c^T B c)/2 = c_average^T B Delta(c).
                pressure.add((0.5 * face.left[i] + 0.5 * face.right[i]) * face.excess_drop[i]);
            }
            let z = self.valences[i];
            let p = self.permeances[i].value();
            if z != 0 && p > 0.0 {
                scale.add((p * face.left[i].max(face.right[i])) * f64::from(z).abs());
                face.positive_carrier |= if z > 0 {
                    face.left[i] > 0.0
                } else {
                    face.right[i] > 0.0
                };
                face.negative_carrier |= if z > 0 {
                    face.right[i] > 0.0
                } else {
                    face.left[i] > 0.0
                };
            }
        }
        pressure.add(problem.trapped_particles_mol[0] / problem.left.solvent_volume_m3());
        pressure.add(-problem.trapped_particles_mol[1] / problem.right.solvent_volume_m3());
        face.osmotic_drop = finite(
            face.rt * pressure.finish()?,
            "original osmotic pressure difference is not representable",
        )?;
        face.current_scale = scale.finish()?;
        Ok(face)
    }

    fn current_at(
        &self,
        face: &Face,
        psi: f64,
        output: &mut [f64; MAX_OSMOTIC_SOLUTES],
    ) -> Result<CurrentPoint> {
        let mut current = Sum::default();
        let mut derivative = Sum::default();
        for (i, &valence) in self.valences.iter().enumerate() {
            let p = self.permeances[i].value();
            if p == 0.0
                || (face.left[i] == 0.0 && face.right[i] == 0.0)
                || (!psi.is_finite() && valence != 0)
            {
                output[i] = 0.0;
                continue;
            }
            let z = f64::from(valence);
            let electric = if valence == 0 { 0.0 } else { z * psi };
            let a = finite(
                face.excess_drop[i] + electric,
                "electrochemical drift is not representable",
            )?;
            let force = if face.left[i] > 0.0 && face.right[i] > 0.0 {
                Some(affinity(
                    face.left[i],
                    face.right[i],
                    face.excess_drop[i],
                    electric,
                )?)
            } else {
                None
            };
            let (flux, slope) = channel_flux(p, a, face.left[i], face.right[i], force)?;
            output[i] = flux;
            current.add(z * flux);
            derivative.add((z * slope) * z);
        }
        Ok(CurrentPoint {
            psi,
            current: current.finish()?,
            derivative: derivative.finish()?,
        })
    }

    fn dissipation(
        &self,
        face: &Face,
        potential: DonnanPotential,
        fluxes: &[f64; MAX_OSMOTIC_SOLUTES],
    ) -> Result<Option<f64>> {
        let mut dissipation = Sum::default();
        let psi = match potential {
            DonnanPotential::Finite {
                dimensionless_left_minus_right,
                ..
            } => Some(dimensionless_left_minus_right),
            DonnanPotential::Indeterminate => Some(0.0),
            _ => None,
        };
        for (i, &z) in self.valences.iter().enumerate() {
            if self.permeances[i].value() == 0.0
                || (psi.is_none() && z != 0)
                || (face.left[i] == 0.0 && face.right[i] == 0.0)
            {
                continue;
            }
            if face.left[i] == 0.0 || face.right[i] == 0.0 {
                return Ok(None);
            }
            let force = affinity(
                face.left[i],
                face.right[i],
                face.excess_drop[i],
                f64::from(z) * psi.unwrap_or(0.0),
            )?;
            let value = finite(fluxes[i] * force, "ionic dissipation is not representable")?;
            if value < 0.0 {
                return Err(ElectrodiffusionError::Invalid(
                    "original electrochemical flux has negative dissipation",
                ));
            }
            dissipation.add(value);
        }
        Ok(Some(finite(
            face.rt * dissipation.finish()?,
            "ionic dissipation is not representable",
        )?))
    }
}

/// Stable original exponential-fit donor streams and their positive dJ/dA.
fn channel_flux(p: f64, a: f64, left: f64, right: f64, force: Option<f64>) -> Result<(f64, f64)> {
    let drift = finite(p * a, "channel drift coefficient is not representable")?;
    let exchange = diffusive_exchange(p, drift)?;
    let forward = finite(
        exchange + drift.max(0.0),
        "forward channel coefficient is not representable",
    )?;
    let reverse = finite(
        exchange + (-drift).max(0.0),
        "reverse channel coefficient is not representable",
    )?;
    let flux = if left == 0.0 {
        -reverse * right
    } else if right == 0.0 {
        forward * left
    } else {
        let force = force.ok_or(ElectrodiffusionError::Invalid(
            "positive channel endpoints need their original affinity",
        ))?;
        // Stream ratios are exp(force). A bounded expm1 preserves detailed
        // balance without subtracting two large, nearly equal streams.
        if force >= 0.0 {
            (forward * left) * -det::expm1(-force)
        } else {
            (reverse * right) * det::expm1(force)
        }
    };
    let magnitude = a.abs();
    let weight = if magnitude < 1e-3 {
        let square = magnitude * magnitude;
        0.5 - magnitude * (1.0 / 6.0 - square * (1.0 / 180.0 - square / 5040.0))
    } else {
        let denominator = -det::expm1(-magnitude);
        det::exp(-magnitude) * (magnitude - denominator) / (denominator * denominator)
    };
    let mean = if a >= 0.0 {
        (1.0 - weight) * left + weight * right
    } else {
        weight * left + (1.0 - weight) * right
    };
    Ok((
        finite(flux, "ionic channel flux is not representable")?,
        finite(p * mean, "ionic current derivative is not representable")?,
    ))
}

// Fluxes and dissipation MUST use the same rounded conjugate affinity. Keeping
// these three terms in one compensated sum also cancels common large forces.
fn affinity(left: f64, right: f64, excess: f64, electric: f64) -> Result<f64> {
    let mut value = Sum::default();
    value.add(log_ratio(left, right));
    value.add(excess);
    value.add(electric);
    Ok(value.finish()?)
}

fn solve_current(
    tolerance: f64,
    mut distance: f64,
    evaluate: &mut impl FnMut(f64) -> Result<CurrentPoint>,
) -> Result<CurrentPoint> {
    let zero = evaluate(0.0)?;
    if zero.current.abs() <= tolerance {
        return Ok(zero);
    }
    let direction = if zero.current > 0.0 { -1.0 } else { 1.0 };
    let mut previous = zero;
    let (mut lo, mut hi) = loop {
        let point = evaluate(direction * distance)?;
        if point.current.abs() <= tolerance {
            return Ok(point);
        }
        if direction * point.current > 0.0 {
            break if direction > 0.0 {
                (previous, point)
            } else {
                (point, previous)
            };
        }
        previous = point;
        distance *= 2.0;
        if !distance.is_finite() {
            return Err(ElectrodiffusionError::NonConvergence);
        }
    };
    let mut bisect_next = false;
    loop {
        let width = hi.psi - lo.psi;
        let best = if lo.current.abs() < hi.current.abs() {
            lo
        } else {
            hi
        };
        let proposed = best.psi - best.current / best.derivative;
        let interpolate =
            !bisect_next && proposed.is_finite() && proposed > lo.psi && proposed < hi.psi;
        let psi = if interpolate {
            proposed
        } else {
            lo.psi + 0.5 * width
        };
        if psi == lo.psi || psi == hi.psi {
            return Err(ElectrodiffusionError::NonConvergence);
        }
        let point = evaluate(psi)?;
        if point.current.abs() <= tolerance {
            return Ok(point);
        }
        if point.current < 0.0 {
            lo = point;
        } else {
            hi = point;
        }
        bisect_next = interpolate && hi.psi - lo.psi > 0.5 * width;
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(ElectrodiffusionError::Cancelled)
    } else {
        Ok(())
    }
}
