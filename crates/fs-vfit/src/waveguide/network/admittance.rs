//! Passive parallel RC memory, through power-preserving pressure/flow duality.
//!
//! Y(s)=G+s*C+sum g*s/(s+p). Each term is a SERIES RC branch in parallel:
//! z'=p*(P-z), U=g*(P-z), H=g*z²/(2p), loss=g*(P-z)².
//! Exchange pressure and flow in the existing RelaxationImpedance midpoint
//! owner: its incident wave is a/Z and its port impedance is 1/Z. This is an
//! algebraic dual, not a new time integrator, delayed thermal force or output EQ.
//! Physical getters keep admittance coefficients and pressure histories distinct
//! from the original impedance's resistance/inertance and flow coordinates.
use crate::impedance::SeriesImpedanceSpec;
use crate::relaxation::{RelaxationFrame, RelaxationImpedance, RelaxationImpedanceSpec, RelaxationTerm, MAX_RELAXATION_BRANCHES};
use crate::waveguide::WaveguideError;

/// One passive series-RC shunt: high-frequency conductance g and relaxation p.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdmittanceTerm {
    /// Positive g [m³/(Pa s)]; branch compliance is g/p [m³/Pa].
    pub conductance_m3_pa_s: f64,
    /// Positive rate [1/s].
    pub rate_per_s: f64,
}
/// Admitted coefficients, retaining the explicit order and all positive terms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxationAdmittanceSpec { dual: RelaxationImpedanceSpec }
impl RelaxationAdmittanceSpec {
    /// Admit nonnegative direct conductance/compliance and up to eight RC arms.
    /// No floor, fit, material identity or omitted-term approximation is implied.
    ///
    /// # Errors
    /// Nonfinite/negative base, nonpositive branch, or exceeded branch budget.
    pub fn new(conductance_m3_pa_s: f64, compliance_m3_pa: f64, terms: &[AdmittanceTerm])
        -> Result<Self, WaveguideError>
    {
        if terms.len()>MAX_RELAXATION_BRANCHES { return Err(WaveguideError("admittance exceeds eight relaxation arms")); }
        let mut dual=[RelaxationTerm::default();MAX_RELAXATION_BRANCHES];
        for (dst,src) in dual.iter_mut().zip(terms) {
            *dst=RelaxationTerm {resistance_pa_s_m3:src.conductance_m3_pa_s,rate_per_s:src.rate_per_s};
        }
        Ok(Self {dual:RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
            resistance_pa_s_m3:conductance_m3_pa_s,inertance_pa_s2_m3:compliance_m3_pa,compliance_m3_pa:None,
        },&dual[..terms.len()])?})
    }
    /// Direct static conductance [m³/(Pa s)].
    #[must_use]
    pub fn conductance_m3_pa_s(&self)->f64 { self.dual.base().resistance_pa_s_m3 }
    /// Direct compliance [m³/Pa], excluding every relaxation arm.
    #[must_use]
    pub fn compliance_m3_pa(&self)->f64 { self.dual.base().inertance_pa_s2_m3 }
    /// Ordered branch coefficients with physical admittance units.
    pub fn terms(&self)->impl ExactSizeIterator<Item=AdmittanceTerm> + '_ {
        self.dual.terms().iter().map(|t|AdmittanceTerm {
            conductance_m3_pa_s:t.resistance_pa_s_m3,rate_per_s:t.rate_per_s,
        })
    }
}
/// Complete immutable trial; hidden coordinates cannot be forged for publication.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdmittanceFrame {
    /// Midpoint common pressure [Pa].
    pub pressure_pa:f64,
    /// Midpoint total flow into the shunt [m³/s].
    pub flow_m3_s:f64,
    /// Reflected pressure wave [Pa].
    pub reflected_pressure_pa:f64,
    /// Actual end-of-step compliance storage [J].
    pub stored_energy_j:f64,
    /// Change in actual storage [J].
    pub storage_change_j:f64,
    /// Irreversible resistor loss [J], not all supplied port work.
    pub dissipated_energy_j:f64,
    /// P*U*dt [J], including possible energy return.
    pub supplied_work_j:f64,
    candidate:RelaxationFrame,
}
impl AdmittanceFrame {
    /// Uncorrected local work defect [J].
    #[must_use]
    pub fn balance_residual_j(&self)->f64 {self.storage_change_j+self.dissipated_energy_j-self.supplied_work_j}
}
/// Bounded inline physical memory; numerical work stays in the existing owner.
#[derive(Clone, Copy, Debug)]
pub struct RelaxationAdmittance { dual:RelaxationImpedance, impedance:f64 }
impl RelaxationAdmittance {
    /// Initially relaxed, zero-energy shunt at a positive pressure/flow impedance.
    ///
    /// # Errors
    /// Invalid port/time or unrepresentable dual coefficients/storage.
    pub fn new(spec:RelaxationAdmittanceSpec, impedance:f64, dt:f64)->Result<Self,WaveguideError> {
        if !impedance.is_finite() || impedance<=0.0 {return Err(WaveguideError("admittance needs a positive finite wave impedance"));}
        Ok(Self {dual:RelaxationImpedance::new(spec.dual,1.0/impedance,dt)?,impedance})
    }
    /// Accepted base compliance pressure [Pa]; zero when no base C is present.
    #[must_use]
    pub fn base_pressure_pa(&self)->f64 {self.dual.base_state().inertive_flow_m3_s}
    /// Accepted internal RC pressures [Pa], not impedance branch flows.
    #[must_use]
    pub fn branch_pressures_pa(&self)->&[f64] {self.dual.branch_flows_m3_s()}
    /// Independently evaluated compliance storage [J].
    #[must_use]
    pub fn stored_energy_j(&self)->f64 {self.dual.stored_energy_j()}
    /// Preview against the current incident wave, without publishing any history.
    ///
    /// # Errors
    /// Nonfinite input, candidate pressure/flow, stored energy or work.
    pub fn preview_step(&self,incident:f64)->Result<AdmittanceFrame,WaveguideError> {
        let candidate=self.dual.preview_step(incident/self.impedance)?;
        let p=candidate.port;
        let result=AdmittanceFrame {pressure_pa:p.flow_m3_s,flow_m3_s:p.pressure_pa,
            reflected_pressure_pa:p.flow_m3_s-incident,stored_energy_j:p.stored_energy_j,
            storage_change_j:p.storage_change_j,dissipated_energy_j:p.dissipated_energy_j,
            supplied_work_j:p.supplied_work_j,candidate};
        if !incident.is_finite() || !result.reflected_pressure_pa.is_finite() {
            return Err(WaveguideError("admittance wave left the finite set"));
        }
        Ok(result)
    }
    pub(crate) fn accept_frame(&mut self,frame:AdmittanceFrame) {self.dual.accept_frame(frame.candidate);}
    /// Advance only a complete accepted trial. Refusal leaves all state intact.
    ///
    /// # Errors
    /// Same admission as [`Self::preview_step`].
    pub fn step(&mut self,incident:f64)->Result<AdmittanceFrame,WaveguideError> {
        let f=self.preview_step(incident)?;self.accept_frame(f);Ok(f)
    }
}

#[cfg(test)]
#[path = "admittance_tests.rs"]
mod tests;
