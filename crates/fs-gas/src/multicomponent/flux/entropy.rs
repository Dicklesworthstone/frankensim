//! Entropy-conjugate inviscid fluxes for the original ideal-mixture calorics.
//!
//! With mathematical entropy eta=-S, beta=1/T and r_i=R/M_i, the entropy
//! variables in the ORIGINAL species-referenced energy chart are
//!   v_i=(h_i,ref/T-s_i)/M_i-beta*u^2/2, v_m=beta*u, v_E=-beta,
//! and the entropy potential is psi=u*sum_i(r_i*rho_i).
//! Here s_i is the molar partial-species entropy, including ideal mixing.
//!
//! For phi_i=beta*u_i,ref-s_i^0+R*ln(T), d(phi_i)/d(beta)=u_i,ref.
//! The divided ORIGINAL internal energy e_i#=Delta(phi_i)/(M_i*Delta beta)
//! therefore replaces the constant-Cv reciprocal logarithmic-temperature mean.
//! The flux uses
//!   F_i=rho_i^log*ubar,
//!   F_m=ubar*sum(F_i)+sum(r_i*rhobar_i)/betabar,
//!   F_E=sum(e_i#*F_i)+(u_left*u_right/2)*sum(F_i)+ubar*p#.
//! Product jump identities give Delta(v).F_ec=Delta(psi) in exact arithmetic.
//! The public receipt retains the actual source-interval and rounded-flux defect.
//!
//! F_es=F_ec-alpha*Delta(q)/2 uses the original endpoint acoustic speeds.
//! Convexity of -S on the admitted common-positive caloric/species domain gives
//! the semidiscrete numerical contribution alpha*Delta(v).Delta(q)/2 >= 0.
//! This is not molecular heating, a positivity timestep theorem, an entropy
//! proof for reconstruction/time integration, or a pressure-equilibrium cure
//! for moving material contacts. Actual time-step admission remains necessary.
//!
//! Calorically perfect multicomponent construction and limitations:
//! Gouasmi, Duraisamy and Murman, https://arxiv.org/abs/1904.00972,
//! Sections 2 and 3. The source-interval divided energy above is derived here
//! for the admitted temperature-dependent, piecewise caloric model.

use super::{
    Cell, Error, Face, FaceFlux, MAX_COMPONENTS, MixtureGas, Primitive, physical, poll, positive,
};
use crate::multicomponent::wall::log_ratio;
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Actual entropy contractions of one oriented inviscid numerical face.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EntropyFaceReport {
    /// Maximum original endpoint |u|+c [m/s]; not a nonlinear positivity bound.
    pub signal_speed: f64,
    /// Delta(psi)-Delta(v).F_ec [W/(m² K)], using the actual original-source
    /// intervals and rounded conservative flux. Signed; never clamped to zero.
    pub entropy_conservative_residual_w_m2_k: f64,
    /// alpha*Delta(v).Delta(q)/2 [W/(m² K)] for the dissipative choice;
    /// zero for the conservative choice. This is numerical entropy production,
    /// not an independent physical heat source. Floating-point sign is retained.
    pub dissipation_entropy_w_m2_k: f64,
    /// Delta(psi)-Delta(v).F_output [W/(m² K)] on the actual returned flux.
    /// Its small difference from the sum above exposes final flux rounding.
    pub entropy_production_w_m2_k: f64,
}

impl MixtureGas {
    /// Evaluate the source-derived entropy-conservative axial mixture flux.
    ///
    /// Output has species_count()+2 entries in the existing species mass,
    /// axial momentum and species-referenced energy order. The physical energy
    /// flux adds each original energy reference times its OWN species flux.
    /// Every active species must be positive on both traces; common absence
    /// is admitted. One-sided zeros refuse without seeding mass or selecting
    /// another flux. Successful face computation allocates no heap storage.
    ///
    /// The entropy identity is semidiscrete and subject to the signed residual
    /// in the report. No positivity, reconstruction or time-step theorem follows.
    /// This central flux alone has no shock dissipation.
    pub fn entropy_conservative_flux(
        &self,
        left: Cell<'_>,
        right: Cell<'_>,
        out: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<EntropyFaceReport, Error> {
        self.entropy_flux::<false>(left, right, out, cancelled)
    }

    /// Evaluate the same source-conjugate flux with scalar acoustic dissipation.
    ///
    /// The subtraction alpha*(q_right-q_left)/2 uses one original |u|+c scale
    /// for ALL species, momentum and referenced energy. It introduces no second
    /// energy/heating source. The measured entropy production is numerical,
    /// with source and final-flux roundoff retained in the returned receipt.
    ///
    /// Common species support is required. The exact-arithmetic semidiscrete
    /// entropy inequality does not guarantee admissible finite-step states or
    /// constant pressure across moving contacts between different caloric laws.
    pub fn entropy_stable_flux(
        &self,
        left: Cell<'_>,
        right: Cell<'_>,
        out: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<EntropyFaceReport, Error> {
        self.entropy_flux::<true>(left, right, out, cancelled)
    }

    fn entropy_flux<const DISSIPATIVE: bool>(
        &self,
        left: Cell<'_>,
        right: Cell<'_>,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EntropyFaceReport, Error> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        let m = n + 2;
        if out.len() != m {
            return Err(Error::Invalid("entropy mixture face output dimension"));
        }
        let pl = self.primitive(left, &mut cancelled)?;
        let pr = self.primitive(right, &mut cancelled)?;
        let mut ql = [0.0; MAX_COMPONENTS];
        let mut qr = [0.0; MAX_COMPONENTS];
        ql[..n].copy_from_slice(left.partial_densities);
        qr[..n].copy_from_slice(right.partial_densities);
        ql[n] = left.momentum;
        qr[n] = right.momentum;
        ql[n + 1] = left.energy;
        qr[n + 1] = right.energy;
        self.interface_entropy_pair::<DISSIPATIVE, true>(
            &ql[..m],
            pl,
            &qr[..m],
            pr,
            out,
            &mut cancelled,
        )
    }

    /// Native first-order owner already admitted these exact source states.
    /// Skip repeated EOS inversion and diagnostic contractions on its hot path;
    /// actual time-step entropy admission remains with the original integrator.
    pub(in crate::multicomponent) fn interface_entropy(
        &self,
        left: &[f64],
        pl: Primitive,
        right: &[f64],
        pr: Primitive,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Face, Error> {
        let report =
            self.interface_entropy_pair::<true, false>(left, pl, right, pr, out, cancelled)?;
        Ok(FaceFlux {
            signal_speed: report.signal_speed,
            used_hll_fallback: false,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn interface_entropy_pair<const DISSIPATIVE: bool, const AUDIT: bool>(
        &self,
        left: &[f64],
        pl: Primitive,
        right: &[f64],
        pr: Primitive,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EntropyFaceReport, Error> {
        poll(cancelled)?;
        let n = self.species_count();
        let m = n + 2;
        if left.len() != m || right.len() != m || out.len() != m {
            return Err(Error::Invalid("entropy mixture face dimensions"));
        }
        let speed = (pl.velocity.abs() + pl.sound_speed).max(pr.velocity.abs() + pr.sound_speed);
        if !positive(speed) || !positive(pl.temperature) || !positive(pr.temperature) {
            return Err(Error::Invalid(
                "entropy mixture face state or acoustic scale",
            ));
        }
        for k in 0..n {
            if (left[k] == 0.0) != (right[k] == 0.0) {
                return Err(Error::Invalid(
                    "entropy mixture flux requires common positive species support",
                ));
            }
        }
        let mut report = EntropyFaceReport {
            signal_speed: speed,
            ..EntropyFaceReport::default()
        };
        if left == right && pl == pr {
            // The authoritative conservative state has already been rounded.
            // Preserve its exact physical flux rather than reconstructing it
            // from source intervals and introducing a uniform-flow defect.
            physical(left, pl, out, n, n + 1);
            check_flux(out)?;
            poll(cancelled)?;
            return Ok(report);
        }

        let half = Dd::from_f64(0.5);
        let tl = Dd::from_f64(pl.temperature);
        let tr = Dd::from_f64(pr.temperature);
        let beta_l = Dd::ONE / tl;
        let beta_r = Dd::ONE / tr;
        let beta_bar = half * beta_l + half * beta_r;
        finite(beta_bar)?;
        // Form the original temperature jump BEFORE reciprocal subtraction.
        // Divide by the larger temperature first to avoid a large intermediate.
        let delta_beta = if pl.temperature == pr.temperature {
            Dd::ZERO
        } else if pr.temperature > pl.temperature {
            -((tr - tl) / tr) / tl
        } else {
            -((tr - tl) / tl) / tr
        };
        if pl.temperature != pr.temperature && finite(delta_beta)? == 0.0 {
            return Err(Error::Invalid("unrepresentable inverse-temperature jump"));
        }
        let ul = Dd::from_f64(pl.velocity);
        let ur = Dd::from_f64(pr.velocity);
        let ubar = half * ul + half * ur;
        let delta_u = ur - ul;
        let kinetic_jump = if AUDIT {
            half * (beta_r * ur * ur - beta_l * ul * ul)
        } else {
            Dd::ZERO
        };
        let sampled = Temperature::new(0.5 * pl.temperature + 0.5 * pr.temperature);
        let initial = Temperature::new(pl.temperature);
        let final_t = Temperature::new(pr.temperature);
        let minimum = Temperature::new(self.domain().minimum_temperature);
        let log_t = Dd::from_f64(log_ratio(pr.temperature, pl.temperature));
        let mut entropies = [Dd::ZERO; MAX_COMPONENTS];
        let mut mass_flux = Dd::ZERO;
        let mut energy_flux = Dd::ZERO;
        let mut gas_density_mean = Dd::ZERO;
        let mut gas_density_jump = Dd::ZERO;

        for (k, model) in self.source().models().iter().enumerate() {
            poll(cancelled)?;
            if left[k] == 0.0 {
                out[k] = 0.0;
                continue;
            }
            let mass = Dd::from_f64(model.molar_mass().value());
            let r = Dd::from_f64(R) / mass;
            let log_density = Dd::from_f64(log_ratio(right[k], left[k]));
            let density_mean = logarithmic_mean(left[k], right[k], log_density)?;
            out[k] = finite(density_mean * ubar)?;
            let actual_flux = Dd::from_f64(out[k]);
            mass_flux = mass_flux + actual_flux;
            let e_sample = Dd::from_f64(
                model
                    .molar_internal_energy_change(minimum, sampled)?
                    .value(),
            );
            let relative = if pl.temperature == pr.temperature {
                Dd::ZERO
            } else {
                Dd::from_f64(
                    model
                        .molar_internal_energy_change(sampled, final_t)?
                        .value(),
                ) / tr
                    - Dd::from_f64(
                        model
                            .molar_internal_energy_change(sampled, initial)?
                            .value(),
                    ) / tl
                    - Dd::from_f64(model.molar_entropy_change(initial, final_t)?.value())
                    + Dd::from_f64(R) * log_t
            };
            let divided_energy = if pl.temperature == pr.temperature {
                e_sample
            } else {
                e_sample + relative / delta_beta
            };
            finite(divided_energy)?;
            energy_flux = energy_flux + divided_energy * (actual_flux / mass);
            let rho_l = Dd::from_f64(left[k]);
            let rho_r = Dd::from_f64(right[k]);
            gas_density_mean = gas_density_mean + r * (half * rho_l + half * rho_r);
            if AUDIT {
                gas_density_jump = gas_density_jump + r * (rho_r - rho_l);
                // Use the actual source interval, not the inferred exact-
                // arithmetic identity divided_energy*delta_beta.
                let caloric_jump = e_sample * delta_beta + relative;
                entropies[k] = caloric_jump / mass + r * log_density - kinetic_jump;
                finite(entropies[k])?;
            }
        }
        let pressure = gas_density_mean / beta_bar;
        out[n] = finite(ubar * mass_flux + pressure)?;
        // This product order avoids squaring a large admitted velocity before
        // it is multiplied by a small mass flux.
        out[n + 1] = finite(energy_flux + (half * ul) * (ur * mass_flux) + ubar * pressure)?;
        check_flux(out)?;

        let mut potential_jump = Dd::ZERO;
        if AUDIT {
            entropies[n] = beta_r * ur - beta_l * ul;
            entropies[n + 1] = -delta_beta;
            potential_jump = ubar * gas_density_jump + delta_u * gas_density_mean;
            let mut contraction = Dd::ZERO;
            for k in 0..m {
                contraction = contraction + entropies[k] * Dd::from_f64(out[k]);
            }
            report.entropy_conservative_residual_w_m2_k = finite(potential_jump - contraction)?;
        }
        if DISSIPATIVE {
            let scale = half * Dd::from_f64(speed);
            let mut contribution = Dd::ZERO;
            for k in 0..m {
                let jump = Dd::from_f64(right[k]) - Dd::from_f64(left[k]);
                if AUDIT {
                    contribution = contribution + entropies[k] * jump;
                }
                out[k] = finite(Dd::from_f64(out[k]) - scale * jump)?;
            }
            if AUDIT {
                report.dissipation_entropy_w_m2_k = finite(scale * contribution)?;
            }
        }
        if AUDIT {
            let mut contraction = Dd::ZERO;
            for k in 0..m {
                contraction = contraction + entropies[k] * Dd::from_f64(out[k]);
            }
            report.entropy_production_w_m2_k = finite(potential_jump - contraction)?;
        }
        poll(cancelled)?;
        Ok(report)
    }
}

// Both arguments are strictly positive original partial densities. The original
// logarithm preserves weak jumps and is also used in the entropy contraction.
fn logarithmic_mean(left: f64, right: f64, logarithm: Dd) -> Result<Dd, Error> {
    if left == right {
        return Ok(Dd::from_f64(left));
    }
    let mean = (Dd::from_f64(right) - Dd::from_f64(left)) / logarithm;
    if !positive(finite(mean)?) {
        return Err(Error::Invalid("unrepresentable logarithmic species mean"));
    }
    Ok(mean)
}

fn finite(value: Dd) -> Result<f64, Error> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(Error::Invalid(
            "non-finite entropy mixture flux or contraction",
        ))
    }
}

fn check_flux(out: &[f64]) -> Result<(), Error> {
    if out.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(Error::Invalid("non-finite entropy mixture flux"))
    }
}

#[cfg(test)]
mod tests;
