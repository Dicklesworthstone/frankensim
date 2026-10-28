//! Scalar, rate-independent magnetic history from an energy-based play model.
//!
//! Each branch stores a polarization z [T], a reversible susceptibility chi
//! [H/m], a pinning field hc [A/m], and a saturation polarization zs [T]. With
//! B = mu_b H + sum(z), the stored energy density and irreversible increment are
//! W = mu_b H^2/2 + sum(z^2/(2 chi)), Q = sum(hc |delta z|).
//! A branch minimizes z^2/(2 chi) - H z + hc |z-z_old| on [-zs,zs].
//! The scalar minimizer is explicit; the B-driven inverse is strictly monotone.
//! This is the uniaxial dry-friction energy construction discussed by Prigozhin
//! et al., IEEE Trans. Magn. 52 (2016), doi:10.1109/TMAG.2016.2599143,
//! <https://arxiv.org/abs/1605.08748>, with a quadratic reversible potential and
//! explicit saturation bounds. It is not their general vector implementation.
//!
//! Supply every field reversal as a separate increment. For a monotone segment,
//! the reported integral of H dB is exact for this piecewise-linear model,
//! including pinning and saturation crossings; endpoint trapezoids are not.
//! Independent branch states retain remanence and minor-loop memory. Parameters
//! are explicit caller fits, not sourced authority inferred from a material name,
//! permeability or a reversible B-H curve. Temperature, stress and model
//! parameters remain fixed. No eddy/excess loss, dynamic hysteresis, vector
//! rotation or permanent-magnet damage law is supplied. Do not add a core-loss
//! map for the same hysteretic heat. This law is opt-in and is deliberately not
//! presented as a path-independent `Storage` function of B alone.

use super::FieldStorageError;
use std::sync::Arc;

/// Isotropic three-dimensional internal polarization and rotating-field memory.
pub mod vector;

/// One polarization contribution; the sum of saturation values is the total
/// saturation polarization. No normalization or hidden branch weight is used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HysteresisBranch {
    /// Reversible slope dz/dH while this branch slides [H/m].
    pub susceptibility_h_m: f64,
    /// Nonnegative dry-friction pinning threshold [A/m].
    pub pinning_field_a_m: f64,
    /// Positive bound |z| <= zs [T].
    pub saturation_polarization_t: f64,
}

#[derive(Debug, PartialEq)]
struct Parameters {
    background_permeability_h_m: f64,
    field_limit_a_m: f64,
    branches: Vec<HysteresisBranch>,
}

/// Immutable, explicitly parameterized uniaxial magnetic material law.
/// The positive background permeability retains an invertible saturated tail.
#[derive(Clone, Debug, PartialEq)]
pub struct ScalarHysteresis {
    parameters: Arc<Parameters>,
}

/// Accepted internal magnetic state. Private fields prevent arbitrary states
/// outside the pinning constraints, or silent replacement of the material law.
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticHysteresisState {
    parameters: Arc<Parameters>,
    field_a_m: f64,
    flux_density_t: f64,
    polarizations_t: Vec<f64>,
    stored_energy_j_m3: f64,
}

impl MagneticHysteresisState {
    /// Accepted H [A/m].
    pub const fn field_a_m(&self) -> f64 {
        self.field_a_m
    }
    /// Accepted B [T], including background and all branch polarizations.
    pub const fn flux_density_t(&self) -> f64 {
        self.flux_density_t
    }
    /// Accepted internal branch polarizations [T].
    pub fn polarizations_t(&self) -> &[f64] {
        &self.polarizations_t
    }
    /// Reversible magnetic energy of this full internal state [J/m^3].
    pub const fn stored_energy_j_m3(&self) -> f64 {
        self.stored_energy_j_m3
    }
}

/// Candidate constitutive increment; the caller commits by retaining `state`.
/// Refusal never changes the previous state or its internal variables.
#[derive(Clone, Debug, PartialEq)]
pub struct HysteresisStep {
    /// Candidate endpoint and complete magnetic memory.
    pub state: MagneticHysteresisState,
    /// Nonnegative pinning heat for this increment [J/m^3].
    pub irreversible_heat_j_m3: f64,
    /// Integral H dB along the admitted monotone segment [J/m^3].
    pub magnetic_work_j_m3: f64,
    /// Stored energy change evaluated without subtracting large energies [J/m^3].
    pub energy_change_j_m3: f64,
    /// Work-conjugate segment mean H; endpoint H for a zero increment [A/m].
    pub mean_field_a_m: f64,
    /// Frozen-previous-state active-set dB/dH [H/m]. At pinning/saturation
    /// kinks the pinned/saturated side is selected; no smoothness is claimed.
    pub incremental_permeability_h_m: f64,
}

fn invalid(what: &'static str) -> FieldStorageError {
    FieldStorageError::Invalid(what)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

impl ScalarHysteresis {
    /// Construct a fixed-condition constitutive model. An empty branch list
    /// admits the linear background alone. At most 1024 branches are admitted;
    /// local inverse work is bounded by 96 safeguarded iterations.
    pub fn new(
        background_permeability_h_m: f64,
        field_limit_a_m: f64,
        branches: Vec<HysteresisBranch>,
    ) -> Result<Self, FieldStorageError> {
        if !positive(background_permeability_h_m)
            || !positive(field_limit_a_m)
            || !(1.0 / background_permeability_h_m).is_finite()
            || branches.len() > 1024
        {
            return Err(invalid(
                "hysteresis needs positive finite background permeability and field domain, with at most 1024 branches",
            ));
        }
        let mut slope = background_permeability_h_m;
        let mut polarization = 0.0;
        let mut energy = 0.5 * (background_permeability_h_m * field_limit_a_m) * field_limit_a_m;
        for branch in &branches {
            let chi = branch.susceptibility_h_m;
            let hc = branch.pinning_field_a_m;
            let zs = branch.saturation_polarization_t;
            if !positive(chi)
                || !positive(zs)
                || !hc.is_finite()
                || hc < 0.0
                || !(zs / chi + hc + field_limit_a_m).is_finite()
            {
                return Err(invalid(
                    "hysteresis branch needs positive susceptibility/saturation and nonnegative finite pinning",
                ));
            }
            slope += chi;
            polarization += zs;
            energy += 0.5 * zs * (zs / chi) + 2.0 * hc * zs;
        }
        if !slope.is_finite()
            || !energy.is_finite()
            || !(background_permeability_h_m * field_limit_a_m + polarization).is_finite()
        {
            return Err(invalid(
                "magnetic history is not representable on its declared field domain",
            ));
        }
        Ok(Self {
            parameters: Arc::new(Parameters {
                background_permeability_h_m,
                field_limit_a_m,
                branches,
            }),
        })
    }

    /// Virgin, unpolarized state at H=B=0. Previously magnetized specimens must
    /// be prepared through their actual field history; no remanence is guessed.
    pub fn initial_state(&self) -> MagneticHysteresisState {
        MagneticHysteresisState {
            parameters: Arc::clone(&self.parameters),
            field_a_m: 0.0,
            flux_density_t: 0.0,
            polarizations_t: vec![0.0; self.parameters.branches.len()],
            stored_energy_j_m3: 0.0,
        }
    }

    fn admit(&self, previous: &MagneticHysteresisState) -> Result<(), FieldStorageError> {
        if !Arc::ptr_eq(&self.parameters, &previous.parameters)
            && self.parameters != previous.parameters
        {
            return Err(invalid(
                "magnetic history must retain the original constitutive parameters",
            ));
        }
        Ok(())
    }

    /// Advance along one monotone H segment. The exact play update retains
    /// saturation and reversal memory without a time-step-dependent loss fit.
    pub fn candidate_field(
        &self,
        previous: &MagneticHysteresisState,
        field_a_m: f64,
    ) -> Result<HysteresisStep, FieldStorageError> {
        self.admit(previous)?;
        if !field_a_m.is_finite() || field_a_m.abs() > self.parameters.field_limit_a_m {
            return Err(invalid(
                "magnetic field is outside the hysteresis model domain",
            ));
        }
        let mu = self.parameters.background_permeability_h_m;
        let (b, slope, _) = self.response(previous, field_a_m);
        let mut polarizations = Vec::with_capacity(self.parameters.branches.len());
        let mut energy = 0.5 * (mu * field_a_m) * field_a_m;
        let mut delta_energy =
            0.5 * mu * (field_a_m - previous.field_a_m) * (field_a_m + previous.field_a_m);
        let mut heat = 0.0;
        for (branch, &old) in self
            .parameters
            .branches
            .iter()
            .zip(&previous.polarizations_t)
        {
            let (next, _) = branch_response(branch, old, field_a_m);
            polarizations.push(next);
            energy += 0.5 * next * (next / branch.susceptibility_h_m);
            delta_energy += 0.5 * (next - old) * ((next + old) / branch.susceptibility_h_m);
            heat += branch.pinning_field_a_m * (next - old).abs();
        }
        let work = delta_energy + heat;
        let db = b - previous.flux_density_t;
        let mean_h = if db == 0.0 {
            if field_a_m != previous.field_a_m {
                return Err(invalid(
                    "magnetic increment is below representable flux resolution",
                ));
            }
            field_a_m
        } else {
            work / db
        };
        if [b, slope, energy, delta_energy, heat, work, mean_h]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid("magnetic history increment is unrepresentable"));
        }
        Ok(HysteresisStep {
            state: MagneticHysteresisState {
                parameters: Arc::clone(&self.parameters),
                field_a_m,
                flux_density_t: b,
                polarizations_t: polarizations,
                stored_energy_j_m3: energy,
            },
            irreversible_heat_j_m3: heat,
            magnetic_work_j_m3: work,
            energy_change_j_m3: delta_energy,
            mean_field_a_m: mean_h,
            incremental_permeability_h_m: slope,
        })
    }

    /// Invert one monotone B segment with a safeguarded active-set Newton solve.
    /// The endpoint reports the recomputed constitutive B, with residual bounded
    /// by 32 ulps of the summed absolute flux contributions. No clipping to the
    /// requested B hides a failed constitutive solve. Only the final candidate
    /// allocates branch state; all inverse iterations reuse the accepted state.
    pub fn candidate_flux_density(
        &self,
        previous: &MagneticHysteresisState,
        flux_density_t: f64,
    ) -> Result<HysteresisStep, FieldStorageError> {
        self.admit(previous)?;
        if !flux_density_t.is_finite() {
            return Err(invalid("magnetic flux density must be finite"));
        }
        if flux_density_t == previous.flux_density_t {
            return self.candidate_field(previous, previous.field_a_m);
        }
        let mut lo = -self.parameters.field_limit_a_m;
        let mut hi = self.parameters.field_limit_a_m;
        let b_lo = self.response(previous, lo).0;
        let b_hi = self.response(previous, hi).0;
        if flux_density_t < b_lo || flux_density_t > b_hi {
            return Err(invalid(
                "magnetic flux requires a field outside the hysteresis domain",
            ));
        }
        let mut h = previous.field_a_m;
        for _ in 0..96 {
            let (b, slope, scale) = self.response(previous, h);
            let residual = b - flux_density_t;
            if residual == 0.0
                || residual.abs() <= 32.0 * f64::EPSILON * (scale + flux_density_t.abs())
            {
                return self.candidate_field(previous, h);
            }
            if residual < 0.0 {
                lo = h;
            } else {
                hi = h;
            }
            let newton = h - residual / slope;
            let next = if newton.is_finite() && newton > lo && newton < hi {
                newton
            } else {
                f64::midpoint(lo, hi)
            };
            if next == h {
                break;
            }
            h = next;
        }
        Err(invalid(
            "magnetic history inverse did not meet its constitutive residual bound",
        ))
    }

    fn response(&self, previous: &MagneticHysteresisState, h: f64) -> (f64, f64, f64) {
        let mut b = self.parameters.background_permeability_h_m * h;
        let mut scale = b.abs();
        let mut slope = self.parameters.background_permeability_h_m;
        for (branch, &old) in self
            .parameters
            .branches
            .iter()
            .zip(&previous.polarizations_t)
        {
            let (z, dz_dh) = branch_response(branch, old, h);
            b += z;
            scale += z.abs();
            slope += dz_dh;
        }
        (b, slope, scale)
    }
}

fn branch_response(branch: &HysteresisBranch, old: f64, h: f64) -> (f64, f64) {
    let chi = branch.susceptibility_h_m;
    let hc = branch.pinning_field_a_m;
    let zs = branch.saturation_polarization_t;
    let reversible = old / chi;
    let driven = if hc == 0.0 {
        h
    } else if h > reversible + hc {
        h - hc
    } else if h < reversible - hc {
        h + hc
    } else {
        return (old, 0.0);
    };
    let saturation_field = zs / chi;
    if driven >= saturation_field {
        (zs, 0.0)
    } else if driven <= -saturation_field {
        (-zs, 0.0)
    } else {
        ((chi * driven).clamp(-zs, zs), chi)
    }
}

/// Uniform fixed-geometry coil/core reduction of the same history law.
/// lambda=N A B, I=l H/N, and all energy/heat densities scale by A l.
/// Signed turns retain winding orientation. Fringing, leakage and copper/eddy
/// losses are separate owners, not hidden adjustments to this core history.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformHystereticCore {
    law: ScalarHysteresis,
    linkage_per_t: f64,
    current_per_field: f64,
    volume_m3: f64,
}

/// Current/flux and energy result for one uniform-core magnetic increment.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticCoreStep {
    /// Constitutive candidate, whose state is retained on acceptance.
    pub material: HysteresisStep,
    /// Total linkage [Wb turns].
    pub flux_linkage_wb: f64,
    /// Endpoint terminal current [A].
    pub current_a: f64,
    /// Integral I d(lambda) divided by delta(lambda) [A].
    pub mean_current_a: f64,
    /// Incremental d(lambda)/dI [H], including the active magnetic branches.
    pub incremental_inductance_h: f64,
    /// Terminal magnetic work [J].
    pub magnetic_work_j: f64,
    /// Pinning heat owed once to the core's thermal owner [J].
    pub core_heat_j: f64,
    /// Complete stored magnetic energy at the endpoint [J].
    pub stored_energy_j: f64,
}

impl UniformHystereticCore {
    /// Bind the constitutive law to a declared uniform path and winding.
    pub fn new(
        law: ScalarHysteresis,
        turns: f64,
        area_m2: f64,
        path_length_m: f64,
    ) -> Result<Self, FieldStorageError> {
        if !turns.is_finite() || turns == 0.0 || !positive(area_m2) || !positive(path_length_m) {
            return Err(invalid(
                "uniform hysteretic core requires finite nonzero turns and positive geometry",
            ));
        }
        let result = Self {
            law,
            linkage_per_t: turns * area_m2,
            current_per_field: path_length_m / turns,
            volume_m3: area_m2 * path_length_m,
        };
        if !result.linkage_per_t.is_finite()
            || result.linkage_per_t == 0.0
            || !result.current_per_field.is_finite()
            || result.current_per_field == 0.0
            || !positive(result.volume_m3)
        {
            return Err(invalid(
                "uniform hysteretic core geometry is unrepresentable",
            ));
        }
        Ok(result)
    }

    /// Virgin core state; preparation histories use the same candidate methods.
    pub fn initial_state(&self) -> MagneticHysteresisState {
        self.law.initial_state()
    }

    /// Actual fixed magnetic volume [m3], also used by a finite core heat cell.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }

    /// Validate retained memory and recover its terminal linkage/current and
    /// stored energy without advancing the magnetic history or generating heat.
    pub fn hold_candidate(
        &self,
        state: &MagneticHysteresisState,
    ) -> Result<HystereticCoreStep, FieldStorageError> {
        self.finish(self.law.candidate_field(state, state.field_a_m())?)
    }

    /// Current-driven uniform core, with explicit irreversible magnetic state.
    pub fn candidate_current(
        &self,
        previous: &MagneticHysteresisState,
        current_a: f64,
    ) -> Result<HystereticCoreStep, FieldStorageError> {
        self.finish(
            self.law
                .candidate_field(previous, current_a / self.current_per_field)?,
        )
    }

    /// Linkage-driven core, suitable for a flux-state circuit's constitutive
    /// residual. Each trial is evaluated from the SAME accepted history until
    /// the complete circuit/thermal step accepts; Newton trials are not history.
    pub fn candidate_linkage(
        &self,
        previous: &MagneticHysteresisState,
        linkage_wb: f64,
    ) -> Result<HystereticCoreStep, FieldStorageError> {
        self.finish(
            self.law
                .candidate_flux_density(previous, linkage_wb / self.linkage_per_t)?,
        )
    }

    fn finish(&self, material: HysteresisStep) -> Result<HystereticCoreStep, FieldStorageError> {
        let result = HystereticCoreStep {
            flux_linkage_wb: self.linkage_per_t * material.state.flux_density_t,
            current_a: self.current_per_field * material.state.field_a_m,
            mean_current_a: self.current_per_field * material.mean_field_a_m,
            incremental_inductance_h: self.linkage_per_t / self.current_per_field
                * material.incremental_permeability_h_m,
            magnetic_work_j: self.volume_m3 * material.magnetic_work_j_m3,
            core_heat_j: self.volume_m3 * material.irreversible_heat_j_m3,
            stored_energy_j: self.volume_m3 * material.state.stored_energy_j_m3,
            material,
        };
        if [
            result.flux_linkage_wb,
            result.current_a,
            result.mean_current_a,
            result.incremental_inductance_h,
            result.magnetic_work_j,
            result.core_heat_j,
            result.stored_energy_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || result.incremental_inductance_h <= 0.0
        {
            return Err(invalid(
                "uniform hysteretic core response is unrepresentable",
            ));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
