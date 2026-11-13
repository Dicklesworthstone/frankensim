//! Pressure-sensitive associated plasticity for a bounded porous solid.
//!
//! Deshpande--Fleck's self-similar ellipse is
//! `phi = sqrt((q^2 + a^2 p^2)/(1 + a^2/9))`, with tension-positive mean
//! stress `p`, von Mises stress `q`, and authored sensitivity `a`. The initial
//! uniaxial yield stress is therefore `yield_stress`. Associated flow retains
//! volumetric plastic strain; hydrostatic compression can leave permanent set.
//! Reference: Deshpande & Fleck (2000), JMPS 48, 1253--1283, Eqs. (2), (3),
//! (11), (13), (21), <https://doi.org/10.1016/S0022-5096(99)00082-4>.
//!
//! This infinitesimal-strain specialization selects constant elastic moduli
//! and linear ENERGETIC hardening. The ellipse is symmetric in mean stress:
//! it is not an asymmetric soil cap, damage law, pore-fluid model, evolving
//! porosity, large-strain densification or an empirical foam calibration.
//! A scalar safeguarded return and its analytic tangent use fixed-size arrays.
//! No finite-difference tangent or tensor Newton iteration is needed.

/// Fixed-condition material sources for the native porous plasticity law.
pub mod source;

pub(super) mod return_map;
use return_map::ReturnParameters as Moduli;

use super::energy::{EnergyTolerances, IncrementEnergy, StoredEnergy};
use crate::elastic::IsotropicElastic;
use crate::tensor::{contract, deviator, trace};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// Authored porous-solid parameters; no density-to-strength fit is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousPlasticity {
    pub elastic: IsotropicElastic,
    /// Initial uniaxial yield stress [Pa].
    pub yield_stress: f64,
    /// Energetic isotropic hardening modulus [Pa], nonnegative.
    pub hardening: f64,
    /// Ellipse parameter `a >= 0`; zero is the pressure-independent J2 limit.
    pub pressure_sensitivity: f64,
}

/// Accepted history. Plastic strain includes its volumetric part.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct PorousState {
    pub plastic_strain: Voigt,
    /// Accumulated work-conjugate plastic strain, not just the Mises strain.
    pub alpha: f64,
}

/// One admitted trial; the input history remains unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousResponse {
    pub stress_pa: Voigt,
    pub tangent_pa: Tangent6,
    pub state: PorousState,
    pub plastic_multiplier: f64,
    pub yield_residual_pa: f64,
}

/// Native response and the shared material work-account representation.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousEnergyIncrement {
    pub strain: Voigt,
    pub response: PorousResponse,
    pub energy: IncrementEnergy,
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("porous plasticity: {what}"),
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

fn hypot(x: f64, y: f64) -> f64 {
    let scale = x.abs().max(y.abs());
    if scale == 0.0 {
        return 0.0;
    }
    scale * det::sqrt((x / scale).powi(2) + (y / scale).powi(2))
}

fn norm(v: &Voigt) -> f64 {
    let scale = v.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    let unit = v.map(|x| x / scale);
    scale * det::sqrt(contract(&unit, &unit))
}

impl PorousPlasticity {
    /// Keep the existing J2 constructor's elastic/yield/hardening ordering.
    pub fn new(
        elastic: IsotropicElastic,
        yield_stress: f64,
        hardening: f64,
        pressure_sensitivity: f64,
    ) -> Result<Self, MaterialError> {
        let law = Self {
            elastic,
            yield_stress,
            hardening,
            pressure_sensitivity,
        };
        law.moduli()?;
        Ok(law)
    }

    fn moduli(&self) -> Result<Moduli, MaterialError> {
        let e = &self.elastic;
        if !positive(e.youngs)
            || !positive(e.strain_limit)
            || !(e.poisson.is_finite() && e.poisson > -1.0 && e.poisson < 0.5)
            || !positive(self.yield_stress)
            || !(self.hardening.is_finite() && self.hardening >= 0.0)
            || !(self.pressure_sensitivity.is_finite() && self.pressure_sensitivity >= 0.0)
        {
            return Err(MaterialError::Parameters {
                what: "porous elasticity, yield, hardening or sensitivity is inadmissible".into(),
            });
        }
        let bulk = e.youngs / (3.0 * (1.0 - 2.0 * e.poisson));
        let shear = e.youngs / (2.0 * (1.0 + e.poisson));
        Moduli::new(
            bulk,
            shear,
            self.yield_stress,
            self.hardening,
            self.pressure_sensitivity,
        )
    }

    fn admit(
        &self,
        strain: &Voigt,
        state: &PorousState,
        m: Moduli,
    ) -> Result<Voigt, MaterialError> {
        let limit = self.elastic.strain_limit.min(0.05);
        if strain.iter().any(|x| !x.is_finite())
            || norm(strain) > limit
            || state.plastic_strain.iter().any(|x| !x.is_finite())
            || !state.alpha.is_finite()
            || state.alpha < 0.0
            || (state.alpha == 0.0 && state.plastic_strain.iter().any(|x| *x != 0.0))
        {
            return Err(invalid(
                "strain or accepted history lies outside its finite small-strain domain",
            ));
        }
        let dev = norm(&deviator(&state.plastic_strain)) * det::sqrt(2.0 / 3.0);
        let volumetric = trace(&state.plastic_strain).abs();
        // Trace roundoff must not be amplified by division by a vanishing a.
        let trace_roundoff = 128.0 * f64::EPSILON * norm(&state.plastic_strain);
        let volume_part = if self.pressure_sensitivity == 0.0 {
            if volumetric > trace_roundoff {
                return Err(invalid("J2 limit requires isochoric history"));
            }
            0.0
        } else {
            (volumetric - trace_roundoff).max(0.0) / self.pressure_sensitivity
        };
        let required = det::sqrt(m.c) * hypot(dev, volume_part);
        if !required.is_finite()
            || required - state.alpha > 128.0 * f64::EPSILON * required.max(state.alpha)
        {
            return Err(invalid(
                "plastic offset exceeds accumulated work-conjugate strain",
            ));
        }
        let elastic = core::array::from_fn(|i| strain[i] - state.plastic_strain[i]);
        if elastic.iter().any(|x| !x.is_finite()) || norm(&elastic) > limit {
            return Err(invalid(
                "elastic strain exceeds the declared small-strain bound",
            ));
        }
        Ok(elastic)
    }

    fn predictor(&self, elastic: &Voigt, m: Moduli) -> (Voigt, f64, f64) {
        let s = deviator(elastic).map(|v| 2.0 * m.shear * v);
        let p = m.bulk * trace(elastic);
        let q = det::sqrt(1.5) * norm(&s);
        (s, p, q)
    }

    /// Checked associated return, including the exact derivative of the discrete update.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &PorousState,
    ) -> Result<PorousResponse, MaterialError> {
        let m = self.moduli()?;
        let elastic = self.admit(strain, state, m)?;
        // Roundoff of the elastic strain (total minus plastic) and its stress at
        // yield must not create flow on an accepted state: a flow below the
        // spacing of alpha cannot advance the history and is refused. The band
        // is the one finite::porous passes to the same return; it is not a
        // physical overstress law.
        let zero = 64.0
            * f64::EPSILON
            * (2.0 * m.shear)
                .max(m.bulk)
                .max(self.yield_stress + self.hardening * state.alpha);
        let returned = m.response(&elastic, state.alpha, zero)?;
        let candidate = if returned.multiplier == 0.0 {
            state.clone()
        } else {
            PorousState {
                plastic_strain: core::array::from_fn(|i| {
                    state.plastic_strain[i] + returned.plastic_increment[i]
                }),
                alpha: returned.alpha,
            }
        };
        self.admit(strain, &candidate, m)?;
        Ok(PorousResponse {
            stress_pa: returned.stress,
            tangent_pa: returned.tangent,
            state: candidate,
            plastic_multiplier: returned.multiplier,
            yield_residual_pa: returned.yield_residual,
        })
    }

    /// Recoverable elastic and explicitly energetic hardening storage [J/m³].
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &PorousState,
    ) -> Result<StoredEnergy, MaterialError> {
        let m = self.moduli()?;
        let e = self.admit(strain, state, m)?;
        let dev = deviator(&e);
        let energy = StoredEnergy {
            elastic_j_m3: 0.5 * m.bulk * trace(&e).powi(2) + m.shear * contract(&dev, &dev),
            hardening_j_m3: 0.5 * self.hardening * state.alpha * state.alpha,
            kinematic_hardening_j_m3: 0.0,
        };
        if !energy.total_j_m3().is_finite() {
            return Err(invalid("stored energy is unrepresentable"));
        }
        Ok(energy)
    }

    /// Audit actual endpoint work without turning integration loss into heat.
    /// The original elastic, yield and energetic hardening coefficients must
    /// be unchanged over the interval. Neither strain nor history is mutated.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &PorousState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<PorousEnergyIncrement, MaterialError> {
        if !positive(tolerances.energy_j_m3) || !positive(tolerances.yield_pa) {
            return Err(invalid(
                "energy and yield allowances must be finite and positive",
            ));
        }
        let m = self.moduli()?;
        let old_e = self.admit(old_strain, old_state, m)?;
        let (_, p, q) = self.predictor(&old_e, m);
        let old_yield = hypot(
            q / det::sqrt(m.c),
            self.pressure_sensitivity / det::sqrt(m.c) * p,
        ) - (self.yield_stress + self.hardening * old_state.alpha);
        if !old_yield.is_finite() || old_yield > tolerances.yield_pa {
            return Err(invalid(
                "old total strain and history are not an admitted committed state",
            ));
        }
        let response = self.response(strain, old_state)?;
        if response.yield_residual_pa > tolerances.yield_pa
            || (response.plastic_multiplier > 0.0
                && response.yield_residual_pa.abs() > tolerances.yield_pa)
        {
            return Err(invalid(
                "endpoint yield residual exceeds the supplied allowance",
            ));
        }
        let new_e = self.admit(strain, &response.state, m)?;
        let (new_s, new_p, _) = self.predictor(&new_e, m);
        let reconstruction: Voigt = core::array::from_fn(|i| {
            new_s[i] + if i < 3 { new_p } else { 0.0 } - response.stress_pa[i]
        });
        if norm(&reconstruction) > tolerances.yield_pa {
            return Err(invalid(
                "committed elastic stress differs from the returned stress",
            ));
        }
        let de: Voigt = core::array::from_fn(|i| new_e[i] - old_e[i]);
        let ds: Voigt = core::array::from_fn(|i| strain[i] - old_strain[i]);
        let da = response.state.alpha - old_state.alpha;
        let dev = deviator(&de);
        let mut energy = IncrementEnergy {
            initial: self.stored_energy(old_strain, old_state)?,
            final_state: self.stored_energy(strain, &response.state)?,
            mechanical_work_j_m3: contract(&response.stress_pa, &ds),
            physical_dissipation_j_m3: self.yield_stress * da,
            numerical_dissipation_j_m3: 0.5 * m.bulk * trace(&de).powi(2)
                + m.shear * contract(&dev, &dev)
                + 0.5 * self.hardening * da * da,
            residual_j_m3: 0.0,
        };
        energy.residual_j_m3 = energy.mechanical_work_j_m3
            - energy.stored_change_j_m3()
            - energy.physical_dissipation_j_m3
            - energy.numerical_dissipation_j_m3;
        if !energy.residual_j_m3.is_finite()
            || energy.residual_j_m3.abs() > tolerances.energy_j_m3
            || !(energy.physical_dissipation_j_m3 >= 0.0
                && energy.numerical_dissipation_j_m3 >= 0.0)
        {
            return Err(invalid(
                "endpoint work account exceeds the supplied allowance",
            ));
        }
        Ok(PorousEnergyIncrement {
            strain: *strain,
            response,
            energy,
        })
    }
}

impl SmallStrainLaw for PorousPlasticity {
    type State = PorousState;
    fn initial_state(&self) -> PorousState {
        PorousState::default()
    }
    fn stress(&self, strain: &Voigt, state: &PorousState) -> Voigt {
        self.response(strain, state)
            .map_or([f64::NAN; 6], |r| r.stress_pa)
    }
    fn tangent(&self, strain: &Voigt, state: &PorousState) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |r| r.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &PorousState) -> PorousState {
        self.response(strain, state)
            .map_or_else(|_| state.clone(), |r| r.state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "normalized Deshpande-Fleck ellipse; bounded small strain",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.porous-plasticity",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain and rotation".into(),
                "associated pressure-sensitive flow".into(),
                "linear energetic hardening; fixed elastic moduli".into(),
                "yield surface symmetric in mean stress".into(),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                self.elastic.strain_limit.min(0.05),
            ),
            vec![
                "no evolving porosity, pore-fluid pressure or large-crush densification".into(),
                "no calibrated material accuracy or asymmetric tensile damage".into(),
            ],
            0.05,
        )
    }
}
