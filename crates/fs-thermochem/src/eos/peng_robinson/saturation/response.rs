//! Homogeneous equilibrium pressure/entropy response inside the coexistence dome.
//!
//! With instantaneous thermal and phase equilibration, `p=p_sat(T)` and
//! `ds=(Cv_eq/T)dT+p_sat' dv`. Therefore `(dp/dv)_s=-T(p_sat')²/Cv_eq`.
//! Cv is the original moving-phase derivative, not frozen-quality calorics.
//! The sound speed additionally assumes common phase velocity and negligible
//! interfacial storage: `c²=-v²(dp/dv)_s/M`. No finite-rate relaxation, slip,
//! attenuation, bubble resonance or nucleation law follows from this limit.
//!
//! Lund and Flåtten (2010), *Equilibrium conditions and sound velocities in
//! two-phase flows*, <https://www.sintef.no/globalassets/project/co2-dynamics/publications/lund_an10.pdf>.

use super::{E, PengRobinsonSaturatedMixtureV1, det};

/// Numerical derivatives of the original pure-fluid equilibrium EOS. The
/// original mixture retains source/phase authority; this snapshot supplies no Cp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PengRobinsonEquilibriumResponseV1 {
    /// Complete molar isochoric heat capacity, including the original ideal Cv.
    pub molar_cv_j_mol_k: f64,
    /// `(dp/dT)_v`, equal to the Clapeyron slope [Pa/K].
    pub dp_dt_v_pa_k: f64,
    /// `(dp/dv)_T`, exactly zero inside the coexistence dome [Pa mol/m³].
    pub dp_dv_t_pa_mol_m3: f64,
    /// `(ds/dT)_v=Cv_eq/T` for molar entropy [J/(mol K²)].
    pub ds_dt_v_j_mol_k2: f64,
    /// `(ds/dv)_T=(dp/dT)_v` for molar entropy [Pa/K].
    pub ds_dv_t_pa_k: f64,
    /// Isentropic pressure-volume derivative [Pa mol/m³].
    pub dp_dv_s_pa_mol_m3: f64,
    /// Homogeneous-equilibrium adiabatic bulk modulus [Pa].
    pub adiabatic_bulk_modulus_pa: f64,
    /// Homogeneous-equilibrium small-signal sound speed [m/s].
    pub sound_speed_m_s: f64,
}

impl PengRobinsonSaturatedMixtureV1 {
    /// Differentiate the original equilibrium at constant volume and entropy.
    /// `ideal_molar_cv_j_mol_k` must come from the same species and temperature.
    /// Exact or numerically unresolved pure-phase endpoints refuse: pressure
    /// response there has distinct one-sided derivatives, not one bidirectional
    /// small-signal sound speed. Both separate phases must be calorically stable.
    pub fn equilibrium_response(
        &self,
        ideal_molar_cv_j_mol_k: f64,
    ) -> Result<PengRobinsonEquilibriumResponseV1, E> {
        let state = self.saturation();
        let volume = self.molar_volume().value();
        if !(0.0 < self.vapor_fraction() && self.vapor_fraction() < 1.0)
            || volume <= state.liquid().molar_volume_m3_per_mol()
            || volume >= state.vapor().molar_volume_m3_per_mol()
        {
            return Err(E::Invalid(
                "bidirectional equilibrium response requires resolved interior phase amount",
            ));
        }
        if !ideal_molar_cv_j_mol_k.is_finite() || ideal_molar_cv_j_mol_k <= 0.0 {
            return Err(E::Invalid("positive original ideal molar Cv required"));
        }
        for root in [state.liquid(), state.vapor()] {
            let cv = ideal_molar_cv_j_mol_k + root.residual_cv_j_mol_k();
            if !cv.is_finite() || cv <= 0.0 {
                return Err(E::Invalid(
                    "both original equilibrium phases require positive complete Cv",
                ));
            }
        }
        let cv = ideal_molar_cv_j_mol_k + self.equilibrium_isochoric_residual_cv_j_mol_k()?;
        let temperature = state.temperature().value();
        let pt = state.saturation_pressure_slope_pa_k();
        let dsdt = cv / temperature;
        let thermal_response = temperature * (pt / cv);
        let dpdv_s = -pt * thermal_response;
        // Avoid v² or p_T² intermediates when each product is representable.
        let bulk = (volume * pt) * thermal_response;
        let sound_squared = bulk / self.mass_density_kg_m3();
        if [cv, pt, dsdt, bulk, sound_squared]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || !dpdv_s.is_finite()
            || dpdv_s >= 0.0
        {
            return Err(E::Unrepresentable(
                "positive equilibrium pressure/entropy response",
            ));
        }
        Ok(PengRobinsonEquilibriumResponseV1 {
            molar_cv_j_mol_k: cv,
            dp_dt_v_pa_k: pt,
            dp_dv_t_pa_mol_m3: 0.0,
            ds_dt_v_j_mol_k2: dsdt,
            ds_dv_t_pa_k: pt,
            dp_dv_s_pa_mol_m3: dpdv_s,
            adiabatic_bulk_modulus_pa: bulk,
            sound_speed_m_s: det::sqrt(sound_squared),
        })
    }
}
