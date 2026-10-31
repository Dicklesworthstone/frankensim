//! Bounded uniaxial thermoviscoelasticity for deposition and thermal histories.
//!
//! Reference `(E_j, tau_j)` pairs have the same meaning as `visco::GeneralizedMaxwell`.
//! Backward Euler is the default; an explicit [`HistoryIntegration::ExactStrainRamp`]
//! selection reuses that model's exact exponential strain-ramp update. Both
//! freeze material coefficients at the proposed endpoint temperature for one
//! interval. Changing temperature therefore still requires time refinement.
//! Temperature changes never reset branch strains.
//! All springs share an explicitly tabulated positive modulus multiplier and
//! relaxation-time multiplier. Free thermal strain is a separate tabulated law.
//!
//! The work audit separates mechanical work, thermal-eigenstrain work,
//! modulus-change work, physical dashpot dissipation, and backward-Euler numerical
//! dissipation. The mechanical update alone is a discrete free-energy balance;
//! `thermal` supplies affine-modulus/expansion entropy and internal energy when explicitly
//! bound to a caloric owner. Mechanical loss is not implicitly deposited as heat.
//! No glass-transition plateau, crystallization, plasticity, damage, finite-strain
//! response, calibration, or three-dimensional constitutive law is inferred.

use std::fmt;

/// Physical dissipation coupled to an existing caloric/heat-transport owner.
pub mod thermal;

/// Source-backed coefficients shared by existing thermal and structural consumers.
pub mod source;

/// Temperature-dependent multipliers and stress-free thermal strain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalKnot {
    /// Absolute temperature, K; knots strictly increase.
    pub temperature_k: f64,
    /// Positive common multiplier of the reference spring moduli.
    pub modulus_scale: f64,
    /// Positive common multiplier of reference relaxation times.
    pub relaxation_scale: f64,
    /// Dimensionless free thermal strain relative to an arbitrary fixed datum.
    /// Each material point subtracts its own value at birth.
    pub free_strain: f64,
}

/// Caller-declared applicability limits; these are not measured error bounds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Limits {
    /// Maximum absolute total AND mechanical strain, at most 0.05.
    pub max_strain: f64,
    /// Maximum total AND mechanical strain-rate magnitude, 1/s.
    pub max_strain_rate_s: f64,
}

/// Constitutive history integration, independent of the declared material data.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum HistoryIntegration {
    /// Endpoint integration with separately accounted algorithmic dissipation.
    #[default]
    BackwardEuler,
    /// Exact Maxwell integration for linear mechanical strain with coefficients
    /// held at the proposed endpoint temperature. Work uses interval-mean stress
    /// and physical dashpot heat is integrated independently; numerical loss is
    /// zero. This is not exact integration of continuously varying temperature.
    ExactStrainRamp,
}

/// A constitutive refusal leaves the accepted point unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Invalid parameters or nonrepresentable arithmetic.
    Invalid(&'static str),
    /// Temperature outside the declared table; no extrapolation.
    Temperature(f64),
    /// Small-strain or strain-rate applicability was exceeded.
    Regime(&'static str),
    /// Actual discrete energy defect, J/m3, exceeds its numerical allowance.
    Energy(f64),
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "thermoviscoelastic refusal: {self:?}")
    }
}
impl std::error::Error for Error {}

/// Immutable generalized-Maxwell law with bounded tabulated thermal response.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermoMaxwell {
    equilibrium_pa: f64,
    branches: Vec<(f64, f64)>,
    knots: Vec<ThermalKnot>,
    limits: Limits,
    history_integration: HistoryIntegration,
}

impl ThermoMaxwell {
    /// Accept the `e_inf` and `terms` of an existing generalized-Maxwell model,
    /// plus explicit thermal data. At most 64 branches and 1024 knots are admitted.
    /// Multipliers and free strain interpolate linearly in temperature. Supply
    /// enough knots to resolve any chosen shift function; none is invented here.
    pub fn new(
        equilibrium_pa: f64,
        branches: Vec<(f64, f64)>,
        knots: Vec<ThermalKnot>,
        limits: Limits,
    ) -> Result<Self, Error> {
        positive(equilibrium_pa)?;
        positive(limits.max_strain)?;
        positive(limits.max_strain_rate_s)?;
        if limits.max_strain > 0.05 || branches.len() > 64 || !(2..=1024).contains(&knots.len()) {
            return Err(Error::Invalid("requires strain <= 0.05, <=64 branches and 2..1024 knots"));
        }
        let mut modulus = equilibrium_pa;
        for &(e, tau) in &branches {
            if !e.is_finite() || e < 0.0 { return Err(Error::Invalid("negative/nonfinite branch modulus")); }
            positive(tau)?;
            modulus = finite(modulus + e)?;
        }
        for knot in &knots {
            positive(knot.temperature_k)?;
            positive(knot.modulus_scale)?;
            positive(knot.relaxation_scale)?;
            finite(knot.free_strain)?;
            positive(finite(modulus * knot.modulus_scale)?)?;
            for &(_, tau) in &branches {
                positive(finite(tau * knot.relaxation_scale)?)?;
            }
        }
        if knots.windows(2).any(|p| p[0].temperature_k >= p[1].temperature_k) {
            return Err(Error::Invalid("temperature knots must strictly increase"));
        }
        Ok(Self { equilibrium_pa, branches, knots, limits,
            history_integration: HistoryIntegration::BackwardEuler })
    }

    /// Select history integration before creating material points. The immutable
    /// law retains every source coefficient and applicability limit; existing
    /// points cannot silently switch their history update.
    #[must_use]
    pub fn with_history_integration(mut self, integration: HistoryIntegration) -> Self {
        self.history_integration = integration;
        self
    }

    #[must_use]
    pub const fn history_integration(&self) -> HistoryIntegration {
        self.history_integration
    }

    /// Number of retained relaxation branches (bounded by construction).
    #[must_use]
    pub fn branch_count(&self) -> usize { self.branches.len() }

    /// Initialize a stress-free material point at its actual deposition strain
    /// and temperature. The point is bound to THIS immutable law for its lifetime.
    pub fn birth(&self, temperature_k: f64, total_strain: f64) -> Result<Point<'_>, Error> {
        let thermal = self.at(temperature_k)?;
        self.check_strain(total_strain)?;
        Ok(Point { law: self, temperature_k, total_strain, birth_strain: total_strain,
            birth_free_strain: thermal.free_strain, q: vec![0.0; self.branches.len()] })
    }

    fn check_strain(&self, strain: f64) -> Result<(), Error> {
        if !strain.is_finite() || strain.abs() > self.limits.max_strain {
            Err(Error::Regime("strain exceeds declared small-strain range"))
        } else { Ok(()) }
    }

    fn at(&self, temperature_k: f64) -> Result<ThermalKnot, Error> {
        if !temperature_k.is_finite() || temperature_k < self.knots[0].temperature_k
            || temperature_k > self.knots[self.knots.len()-1].temperature_k {
            return Err(Error::Temperature(temperature_k));
        }
        let upper = self.knots.partition_point(|k| k.temperature_k < temperature_k);
        if upper == 0 { return Ok(self.knots[0]); }
        let a = self.knots[upper-1];
        let b = self.knots[upper];
        let weight = (temperature_k-a.temperature_k)/(b.temperature_k-a.temperature_k);
        let blend = |x: f64, y: f64| (1.0-weight)*x + weight*y;
        Ok(ThermalKnot { temperature_k,
            modulus_scale: blend(a.modulus_scale, b.modulus_scale),
            relaxation_scale: blend(a.relaxation_scale, b.relaxation_scale),
            free_strain: blend(a.free_strain, b.free_strain) })
    }
}

/// Accepted history. Private state cannot be reassigned to another material law.
#[derive(Debug, Clone, PartialEq)]
pub struct Point<'a> {
    law: &'a ThermoMaxwell,
    temperature_k: f64,
    total_strain: f64,
    birth_strain: f64,
    birth_free_strain: f64,
    // Elastic strains of Maxwell springs, as in GeneralizedMaxwell/PronyState.
    q: Vec<f64>,
}

/// Affine endpoint stress for an implicit structural solve: sigma=C*epsilon+b.
#[derive(Debug, Clone, Copy)]
pub struct AffineStress {
    /// Exact derivative with respect to the new TOTAL strain, Pa.
    pub tangent_pa: f64,
    /// History, birth reference and thermal eigenstrain contribution, Pa.
    pub offset_pa: f64,
}

/// Work/energy densities for the actual accepted interval, all in J/m3.
#[derive(Debug, Clone, Copy, Default)]
pub struct Energy {
    pub mechanical_work: f64,
    pub thermal_strain_work: f64,
    pub modulus_work: f64,
    pub stored_change: f64,
    /// Nonnegative physical dashpot dissipation of the selected history update.
    pub dissipated: f64,
    /// Nonnegative algorithmic dissipation; zero for exact strain ramps.
    pub numerical_dissipation: f64,
    /// Storage + both dissipations - all three work contributions.
    pub residual: f64,
}

/// Proposed next point and independently evaluated interval balance.
#[derive(Debug, Clone)]
pub struct Step<'a> {
    pub point: Point<'a>,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub energy: Energy,
}

impl<'a> Point<'a> {
    #[must_use]
    pub fn temperature_k(&self) -> f64 { self.temperature_k }
    #[must_use]
    pub fn total_strain(&self) -> f64 { self.total_strain }
    #[must_use]
    pub fn birth_strain(&self) -> f64 { self.birth_strain }
    #[must_use]
    pub fn branch_strains(&self) -> &[f64] { &self.q }
    #[must_use]
    pub fn law(&self) -> &'a ThermoMaxwell { self.law }

    fn mechanical_strain(&self) -> f64 {
        self.total_strain-self.birth_strain
            -(self.law.at(self.temperature_k).expect("accepted temperature").free_strain-self.birth_free_strain)
    }
    fn unscaled_storage(&self) -> f64 {
        let e = self.mechanical_strain();
        0.5*self.law.equilibrium_pa*e*e + self.law.branches.iter().zip(&self.q)
            .map(|(&(modulus, _), &q)| 0.5*modulus*q*q).sum::<f64>()
    }
    #[must_use]
    /// Spring potential (Helmholtz free energy in the admitted caloric model).
    pub fn stored_energy_j_m3(&self) -> f64 {
        self.law.at(self.temperature_k).expect("accepted temperature").modulus_scale*self.unscaled_storage()
    }
    #[must_use]
    pub fn stress_pa(&self) -> f64 {
        self.law.at(self.temperature_k).expect("accepted temperature").modulus_scale
            * (self.law.equilibrium_pa*self.mechanical_strain() + self.law.branches.iter().zip(&self.q)
                .map(|(&(e, _), &q)| e*q).sum::<f64>())
    }

    /// Exact affine stress for a positive-duration step. Does not age state.
    /// Applicability of the unknown strain/rate is checked on the final proposal.
    pub fn linearize(&self, temperature_k: f64, duration_s: f64) -> Result<AffineStress, Error> {
        positive(duration_s)?;
        let thermal = self.law.at(temperature_k)?;
        let old_e = self.mechanical_strain();
        let mut tangent = self.law.equilibrium_pa;
        let mut history = 0.0;
        for (&(e, tau), &q) in self.law.branches.iter().zip(&self.q) {
            match self.law.history_integration {
                HistoryIntegration::BackwardEuler => {
                    let (a, _) = factors(tau*thermal.relaxation_scale, duration_s);
                    tangent += e*a;
                    history += e*a*(q-old_e);
                }
                HistoryIntegration::ExactStrainRamp => {
                    let ramp = crate::visco::PronyRamp::new(duration_s/(tau*thermal.relaxation_scale));
                    tangent += e*ramp.apply(0.0, 1.0).0;
                    history += e*ramp.apply(q, -old_e).0;
                }
            }
        }
        tangent = finite(tangent*thermal.modulus_scale)?;
        let offset = finite(thermal.modulus_scale*history
            - tangent*(self.birth_strain+thermal.free_strain-self.birth_free_strain))?;
        positive(tangent)?;
        Ok(AffineStress { tangent_pa: tangent, offset_pa: offset })
    }

    /// Propose the selected history update without mutating accepted state.
    /// Exact strain ramps use mean stress for total/eigenstrain work and
    /// independently integrated physical heat. The temperature-dependent
    /// coefficient change is still held at the proposed endpoint temperature.
    pub fn advance(&self, temperature_k: f64, total_strain: f64, duration_s: f64,
        energy_tolerance_j_m3: f64) -> Result<Step<'a>, Error> {
        positive(energy_tolerance_j_m3)?;
        positive(duration_s)?;
        let backward_affine = if self.law.history_integration == HistoryIntegration::BackwardEuler {
            Some(self.linearize(temperature_k, duration_s)?)
        } else { None };
        self.law.check_strain(total_strain)?;
        let old = self.law.at(self.temperature_k)?;
        let thermal = self.law.at(temperature_k)?;
        let e_old = self.mechanical_strain();
        let d_thermal = thermal.free_strain-old.free_strain;
        let d_total = total_strain-self.total_strain;
        let de = d_total-d_thermal;
        let e_new = total_strain-self.birth_strain-(thermal.free_strain-self.birth_free_strain);
        self.law.check_strain(e_new)?;
        if !de.is_finite() || d_total.abs()/duration_s > self.law.limits.max_strain_rate_s
            || de.abs()/duration_s > self.law.limits.max_strain_rate_s {
            return Err(Error::Regime("strain rate exceeds declared material range"));
        }
        let mut next = self.clone();
        next.temperature_k = temperature_k;
        next.total_strain = total_strain;
        let mut dissipated = 0.0;
        let mut numerical = 0.0;
        let mut mean_stress = self.law.equilibrium_pa*(e_old+0.5*(e_new-e_old));
        let mut exact_tangent = self.law.equilibrium_pa;
        match self.law.history_integration {
            HistoryIntegration::BackwardEuler => {
                numerical = 0.5*self.law.equilibrium_pa*de*de;
                for (j, &(e, tau)) in self.law.branches.iter().enumerate() {
                    let (a, b) = factors(tau*thermal.relaxation_scale, duration_s);
                    let drive = self.q[j] + (e_new-e_old);
                    next.q[j] = finite(a*drive)?;
                    dissipated += (e*a*drive)*(b*drive);
                    let dq = next.q[j]-self.q[j];
                    numerical += 0.5*e*dq*dq;
                }
            }
            HistoryIntegration::ExactStrainRamp => {
                for (j, &(e, tau)) in self.law.branches.iter().enumerate() {
                    let ramp = crate::visco::PronyRamp::new(duration_s/(tau*thermal.relaxation_scale));
                    let (q, mean, heat) = ramp.apply(self.q[j], e_new-e_old);
                    next.q[j] = finite(q)?;
                    exact_tangent += e*ramp.apply(0.0, 1.0).0;
                    mean_stress += e*mean;
                    dissipated += e*heat;
                }
            }
        }
        let stress = finite(next.stress_pa())?;
        let work_stress = match self.law.history_integration {
            HistoryIntegration::BackwardEuler => stress,
            HistoryIntegration::ExactStrainRamp => finite(thermal.modulus_scale*mean_stress)?,
        };
        let tangent_pa = match backward_affine {
            Some(affine) => affine.tangent_pa,
            None => finite(thermal.modulus_scale*exact_tangent)?,
        };
        positive(tangent_pa)?;
        let mut energy = Energy {
            mechanical_work: finite(work_stress*d_total)?,
            thermal_strain_work: finite(-work_stress*d_thermal)?,
            modulus_work: finite((thermal.modulus_scale-old.modulus_scale)*self.unscaled_storage())?,
            stored_change: finite(next.stored_energy_j_m3()-self.stored_energy_j_m3())?,
            dissipated: finite(thermal.modulus_scale*dissipated)?,
            numerical_dissipation: finite(thermal.modulus_scale*numerical)?,
            residual: 0.0,
        };
        energy.residual = finite(energy.stored_change+energy.dissipated+energy.numerical_dissipation
            -energy.mechanical_work-energy.thermal_strain_work-energy.modulus_work)?;
        if energy.residual.abs() > energy_tolerance_j_m3 {
            return Err(Error::Energy(energy.residual));
        }
        Ok(Step { point: next, stress_pa: stress, tangent_pa, energy })
    }
}

// Complementary factors calculated without overflowing tau+dt or losing dt/tau
// in 1-a. Both underflow limits are the corresponding elastic/relaxed limits.
fn factors(tau: f64, dt: f64) -> (f64, f64) {
    if dt <= tau { let r = dt/tau; (1.0/(1.0+r), r/(1.0+r)) }
    else { let r = tau/dt; (r/(1.0+r), 1.0/(1.0+r)) }
}
fn positive(x: f64) -> Result<(), Error> {
    if x.is_finite() && x > 0.0 { Ok(()) } else { Err(Error::Invalid("positive finite quantity required")) }
}
fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() { Ok(x) } else { Err(Error::Invalid("nonfinite constitutive arithmetic")) }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn law(branches: Vec<(f64, f64)>) -> ThermoMaxwell {
        ThermoMaxwell::new(1e8, branches, vec![
            ThermalKnot { temperature_k: 300.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
            ThermalKnot { temperature_k: 500.0, modulus_scale: 0.5, relaxation_scale: 0.1, free_strain: 0.002 },
        ], Limits { max_strain: 0.03, max_strain_rate_s: 10.0 }).unwrap()
    }
    fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol, "{a:e} != {b:e}"); }
    #[test]
    fn free_contraction_and_constrained_cooling_have_distinct_work() {
        let law = law(vec![]); let born = law.birth(500.0, 0.0).unwrap();
        let free = born.advance(300.0, -0.002, 1.0, 1e-8).unwrap();
        close(free.stress_pa, 0.0, 1e-8);
        let held = born.advance(300.0, 0.0, 1.0, 1e-8).unwrap();
        close(held.stress_pa, 2e5, 1e-8);
        close(held.energy.mechanical_work, 0.0, 1e-8);
        close(held.energy.thermal_strain_work, 400.0, 1e-8);
        close(held.energy.stored_change, 200.0, 1e-8);
        close(held.energy.numerical_dissipation, 200.0, 1e-8);
    }
    #[test]
    fn isothermal_maxwell_relaxation_matches_discrete_reference_and_refines() {
        let law = law(vec![(3e8, 2.0)]);
        let loaded = law.birth(300.0, 0.0).unwrap().advance(300.0, 0.01, 0.01, 1e-7).unwrap().point;
        let q0 = loaded.q[0]; let mut errors = Vec::new();
        for n in [10, 20, 40] {
            let mut point = loaded.clone(); let dt = 2.0/f64::from(n);
            for _ in 0..n { point = point.advance(300.0, 0.01, dt, 1e-7).unwrap().point; }
            close(point.q[0], q0*(1.0+dt/2.0).powi(-n), 1e-14);
            errors.push((point.q[0]-q0*(-1.0_f64).exp()).abs());
        }
        assert!(errors[0]/errors[1] > 1.9 && errors[1]/errors[2] > 1.9);
    }
    #[test]
    fn temperature_changes_scale_stiffness_without_erasing_branch_memory() {
        let law = law(vec![(3e8, 2.0)]);
        let loaded = law.birth(300.0, 0.0).unwrap().advance(300.0, 0.01, 0.01, 1e-7).unwrap().point;
        // Match free expansion so mechanical strain stays fixed during heating.
        let hot = loaded.advance(500.0, 0.012, 0.2, 1e-7).unwrap();
        close(hot.point.q[0], 0.5*loaded.q[0], 1e-14);
        assert!(hot.energy.modulus_work < 0.0 && hot.energy.dissipated > 0.0);
        let cold = loaded.advance(300.0, 0.01, 0.2, 1e-7).unwrap();
        assert!(cold.point.q[0] > hot.point.q[0]);
        assert_eq!(loaded.temperature_k(), 300.0);
    }
    #[test]
    fn affine_tangent_matches_actual_endpoint_and_finite_difference() {
        let law = law(vec![(3e8, 2.0), (2e8, 10.0)]);
        let old = law.birth(450.0, 0.001).unwrap().advance(400.0, 0.002, 0.3, 1e-7).unwrap().point;
        let affine = old.linearize(350.0, 0.4).unwrap();
        let a = old.advance(350.0, 0.003, 0.4, 1e-7).unwrap();
        let b = old.advance(350.0, 0.003+1e-7, 0.4, 1e-7).unwrap();
        close(a.stress_pa, affine.tangent_pa*0.003+affine.offset_pa, 1e-6);
        close((b.stress_pa-a.stress_pa)/1e-7/affine.tangent_pa, 1.0, 1e-8);
    }
    #[test]
    fn cyclic_cooling_reheating_closes_all_energy_accounts() {
        let law = law(vec![(3e8, 0.1), (2e8, 10.0)]);
        let mut point = law.birth(500.0, 0.0).unwrap();
        for i in 0..200 {
            let t = if i < 100 { 500.0-2.0*f64::from(i) } else { 300.0+2.0*f64::from(i-100) };
            let step = point.advance(t, 0.0, 0.1, 1e-7).unwrap();
            assert!(step.energy.dissipated >= 0.0 && step.energy.numerical_dissipation >= 0.0);
            assert!(step.energy.residual.abs() < 1e-7); point = step.point;
        }
        assert!(point.branch_strains().iter().any(|q| q.abs() > 1e-7));
    }
    #[test]
    fn domain_and_rate_refusals_preserve_accepted_history() {
        let law = law(vec![(3e8, 1.0)]); let point = law.birth(400.0, 0.0).unwrap();
        let original = point.clone();
        assert!(point.advance(299.0, 0.0, 1.0, 1e-7).is_err());
        assert!(point.advance(400.0, 0.04, 1.0, 1e-7).is_err());
        assert!(point.advance(400.0, 0.01, 1e-8, 1e-7).is_err());
        assert!(point.advance(400.0, 0.01, 0.0, 1e-7).is_err());
        assert_eq!(point, original);
    }
    #[test]
    fn extreme_positive_time_scales_have_finite_limit_updates() {
        let law = law(vec![(3e8, 1e-200), (2e8, 1e200)]);
        let step = law.birth(300.0, 0.0).unwrap().advance(300.0, 0.01, 1.0, 1e-7).unwrap();
        close(step.point.q[0], 0.0, 1e-100); close(step.point.q[1], 0.01, 1e-14);
        assert!(step.energy.residual.abs() < 1e-7);
    }

    #[test]
    fn exact_history_matches_the_original_prony_owner_on_reversing_ramps() {
        let law = law(vec![(3e8, 0.1), (2e8, 10.0)])
            .with_history_integration(HistoryIntegration::ExactStrainRamp);
        let reference = crate::visco::GeneralizedMaxwell::new(1e8, law.branches.clone()).unwrap();
        let mut native = reference.state();
        let mut point = law.birth(300.0, 0.0).unwrap();
        for (strain, dt) in [(0.01, 0.01), (-0.005, 0.2), (0.002, 20.0), (0.002, 1e200)] {
            let (old_work, old_heat) = (native.work, native.dissipated);
            let stress = reference.step(&mut native, strain, dt);
            let step = point.advance(300.0, strain, dt, 1e-7).unwrap();
            close(step.stress_pa, stress, 1e-8);
            for (&actual, &expected) in step.point.q.iter().zip(&native.q) {
                close(actual, expected, 1e-16);
            }
            close(step.energy.mechanical_work, native.work-old_work, 1e-8);
            close(step.energy.dissipated, native.dissipated-old_heat, 1e-8);
            assert_eq!(step.energy.numerical_dissipation, 0.0);
            assert!(step.energy.dissipated >= 0.0);
            point = step.point;
        }
    }

    #[test]
    fn exact_hold_retains_the_original_relaxation_time_under_partition() {
        let law = law(vec![(3e8, 0.001)])
            .with_history_integration(HistoryIntegration::ExactStrainRamp);
        let loaded = law.birth(300.0, 0.0).unwrap()
            .advance(300.0, 0.001, 0.001, 1e-7).unwrap().point;
        let expected = loaded.q[0]*(-10.0_f64).exp();
        let expected_heat = 0.5*3e8*(loaded.q[0]*loaded.q[0]-expected*expected);
        for count in [1, 4, 16, 64] {
            let mut point = loaded.clone();
            let mut heat = 0.0;
            for _ in 0..count {
                let step = point.advance(300.0, 0.001, 0.01/f64::from(count), 1e-7).unwrap();
                assert!(step.point.q[0] >= 0.0 && step.point.q[0] <= point.q[0]);
                assert_eq!(step.energy.numerical_dissipation, 0.0);
                heat += step.energy.dissipated;
                point = step.point;
            }
            close(point.q[0], expected, 1e-18);
            close(heat, expected_heat, 1e-10);
        }
    }

    #[test]
    fn exact_thermal_ramp_has_mean_work_and_the_matching_endpoint_tangent() {
        let spring = law(vec![]).with_history_integration(HistoryIntegration::ExactStrainRamp);
        let held = spring.birth(500.0, 0.0).unwrap().advance(300.0, 0.0, 1.0, 1e-8).unwrap();
        close(held.energy.thermal_strain_work, 200.0, 1e-10);
        close(held.energy.stored_change, 200.0, 1e-10);
        assert_eq!(held.energy.numerical_dissipation, 0.0);

        let law = law(vec![(3e8, 2.0), (2e8, 10.0)])
            .with_history_integration(HistoryIntegration::ExactStrainRamp);
        let old = law.birth(450.0, 0.001).unwrap().advance(400.0, 0.002, 0.3, 1e-7).unwrap().point;
        let affine = old.linearize(350.0, 0.4).unwrap();
        let a = old.advance(350.0, 0.003, 0.4, 1e-7).unwrap();
        let b = old.advance(350.0, 0.003+1e-7, 0.4, 1e-7).unwrap();
        close(a.stress_pa, affine.tangent_pa*0.003+affine.offset_pa, 1e-6);
        close(a.tangent_pa, affine.tangent_pa, 1e-6);
        close((b.stress_pa-a.stress_pa)/1e-7/affine.tangent_pa, 1.0, 1e-8);
        assert!(a.energy.modulus_work > 0.0 && a.energy.dissipated > 0.0);
        assert_eq!(a.energy.numerical_dissipation, 0.0);
    }

    #[test]
    fn exact_history_refusals_leave_accepted_state_unchanged() {
        let law = law(vec![(3e8, 1.0)])
            .with_history_integration(HistoryIntegration::ExactStrainRamp);
        let point = law.birth(400.0, 0.0).unwrap();
        let original = point.clone();
        for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(point.advance(400.0, 0.001, dt, 1e-7).is_err());
        }
        assert!(point.advance(299.0, 0.0, 1.0, 1e-7).is_err());
        assert!(point.advance(400.0, 0.04, 1.0, 1e-7).is_err());
        assert!(point.advance(400.0, 0.01, 1e-8, 1e-7).is_err());
        assert_eq!(point, original);
    }
}
