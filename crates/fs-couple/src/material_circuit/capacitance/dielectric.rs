//! Source-backed nonlinear dielectric charge laws on the ordinary native circuit.
//!
//! Fixed area A and gap g map every original source knot exactly to q=A*D and
//! V=g*E. The existing native V(q) curve owns its integral, discrete gradient
//! and charge history. This module retains the original material and geometry;
//! it creates no fitted table, second energy store or alternate circuit solver.

use fs_material::field_storage::dielectric::reversible::{
    ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY, ReversibleDielectricMaterial,
};
use fs_phs::circuit::capacitance::{CapacitorChargeCurve, CapacitorCurve, MAX_CAPACITOR_KNOTS};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};

use super::{finite, invalid};
use crate::material_circuit::{
    MAX_MATERIAL_CIRCUIT_BRANCHES, MaterialCircuit, MaterialCircuitError, MaterialResistor,
    contacts::SeriesElectricalContact, poll, storage::MaterialStorageCandidate,
};

/// One immutable constitutive source and actual uniform-gap geometry.
/// The prepared native curve and material path both share their data on clone.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialDielectricCapacitor {
    branch: usize,
    material: ReversibleDielectricMaterial,
    area_m2: f64,
    gap_m: f64,
    curve: CapacitorChargeCurve,
}

impl MaterialDielectricCapacitor {
    /// Preserve every source knot and mirror its isotropic zero-remanence law
    /// to both charge polarities. No knot is fitted, dropped or extrapolated.
    /// The native curve's existing knot budget bounds the doubled source path.
    ///
    /// A closed ordinary circuit requires a field-only source potential.
    /// Actual temperature or other external response dependence refuses:
    /// freezing such a free-energy snapshot does not account for its missing
    /// reversible thermal/mechanical exchange. Fixed applicability bounds and
    /// their original query context remain attached to the material.
    pub fn new(
        branch: usize,
        material: ReversibleDielectricMaterial,
        area_m2: f64,
        gap_m: f64,
    ) -> Result<Self, MaterialCircuitError> {
        if branch >= MAX_MATERIAL_CIRCUIT_BRANCHES
            || [area_m2, gap_m].iter().any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "bounded dielectric branch and positive finite area/gap required",
            ));
        }
        let property = material
            .material()
            .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
            .ok_or_else(|| invalid("nonlinear capacitor lost its original dielectric source"))?;
        let axes = property.source_axes();
        if axes.len() != 1
            || axes[0].0 != material.source_path().axis
            || axes[0].1 != ELECTRIC_FIELD_DIMS
        {
            return Err(invalid(
                "closed nonlinear capacitor requires a field-only source potential",
            ));
        }
        let path = &material.source_path().knots;
        let count = path
            .len()
            .checked_mul(2)
            .and_then(|n| n.checked_sub(1))
            .filter(|&n| n <= MAX_CAPACITOR_KNOTS)
            .ok_or_else(|| {
                invalid("original dielectric path exceeds native signed capacitor knot capacity")
            })?;
        let mut knots = Vec::with_capacity(count);
        for &(field, ref answer) in path.iter().skip(1).rev() {
            knots.push((
                -finite(area_m2 * answer.evidence.value.value)?,
                -finite(gap_m * field)?,
            ));
        }
        for &(field, ref answer) in path {
            knots.push((
                finite(area_m2 * answer.evidence.value.value)?,
                finite(gap_m * field)?,
            ));
        }
        let curve = CapacitorChargeCurve::try_new(knots)?;
        Ok(Self {
            branch,
            material,
            area_m2,
            gap_m,
            curve,
        })
    }

    pub const fn branch(&self) -> usize {
        self.branch
    }
    pub fn material(&self) -> &ReversibleDielectricMaterial {
        &self.material
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    /// The exact original native charge law used by every circuit proposal.
    pub fn curve(&self) -> &CapacitorChargeCurve {
        &self.curve
    }
    pub fn maximum_charge_c(&self) -> f64 {
        self.curve.charge_range_c().1
    }
    pub fn maximum_voltage_v(&self) -> f64 {
        self.curve.knots().last().unwrap().1
    }
    /// Signed terminal voltage from the retained native charge law.
    pub fn voltage_v(&self, charge_c: f64) -> Result<f64, MaterialCircuitError> {
        Ok(self.curve.voltage_v(charge_c)?)
    }
    /// Signed initial charge at an admitted voltage on that same law.
    pub fn charge_c(&self, voltage_v: f64) -> Result<f64, MaterialCircuitError> {
        Ok(self.curve.charge_c(voltage_v)?)
    }
    /// Native integral V dq = A*g*integral E dD, without a second energy owner.
    pub fn stored_energy_j(&self, charge_c: f64) -> Result<f64, MaterialCircuitError> {
        Ok(self.curve.stored_energy_j(charge_c)?)
    }
    fn binding(&self) -> CapacitorCurve {
        CapacitorCurve {
            branch: self.branch,
            curve: self.curve.clone(),
        }
    }
}

fn canonical_sources(
    sources: &mut [MaterialDielectricCapacitor],
) -> Result<(), MaterialCircuitError> {
    if sources.len() > MAX_MATERIAL_CIRCUIT_BRANCHES {
        return Err(invalid("dielectric bindings exceed native branch capacity"));
    }
    sources.sort_by_key(MaterialDielectricCapacitor::branch);
    if sources
        .windows(2)
        .any(|pair| pair[0].branch == pair[1].branch)
    {
        return Err(invalid(
            "each dielectric source must own a distinct original capacitor branch",
        ));
    }
    Ok(())
}

impl MaterialCircuit {
    /// Build an ordinary circuit directly from complete nonlinear dielectric
    /// sources and fixed geometries, alongside other authored charge curves,
    /// bulk/contact loss owners, mutual inductance and original source ports.
    /// The selected graph coefficients are lowered to the source's zero-field
    /// differential capacitance; the native curve owns its actual full response.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_dielectric_capacitors(
        mut graph: CircuitGraph,
        resistors: Vec<MaterialResistor>,
        contacts: Vec<SeriesElectricalContact>,
        mutual: Vec<MutualInductance>,
        mut curves: Vec<CapacitorCurve>,
        mut dielectrics: Vec<MaterialDielectricCapacitor>,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        poll(&mut cancelled)?;
        canonical_sources(&mut dielectrics)?;
        if curves
            .len()
            .checked_add(dielectrics.len())
            .is_none_or(|n| n > graph.branches.len() || n > MAX_MATERIAL_CIRCUIT_BRANCHES)
        {
            return Err(invalid("capacitor laws exceed native branch capacity"));
        }
        for source in &dielectrics {
            poll(&mut cancelled)?;
            if curves.iter().any(|curve| curve.branch == source.branch) {
                return Err(invalid(
                    "an authored curve cannot also own a dielectric source branch",
                ));
            }
            let branch = graph.branches.get_mut(source.branch).ok_or_else(|| {
                invalid("dielectric source must name an original capacitor branch")
            })?;
            if !matches!(branch.2, Branch::Capacitor { .. }) {
                return Err(invalid(
                    "dielectric source cannot replace a non-capacitor branch",
                ));
            }
            branch.2 = Branch::Capacitor {
                farads: source.curve.differential_capacitance_f(0.0)?,
            };
            curves.push(source.binding());
        }
        let mut circuit = Self::new_with_capacitor_curves(
            graph, resistors, contacts, mutual, curves, fluxes_wb, charges_c, inputs,
        )?;
        circuit.dielectric_capacitors = dielectrics.into();
        poll(&mut cancelled)?;
        Ok(circuit)
    }

    /// Complete immutable source/geometry records, retained through ordinary
    /// resistor, contact, conductivity and thermal candidate rebuilds.
    pub fn dielectric_capacitors(&self) -> &[MaterialDielectricCapacitor] {
        &self.dielectric_capacitors
    }

    /// Replace every retained dielectric source at fixed original geometry,
    /// charge, flux and topology. Other authored curves/material records and
    /// mutual fields survive unchanged. Source substitution is explicit; the
    /// candidate retains each proposed card and reports signed parameter work.
    ///
    /// The original union-of-knots integral of V_new-V_old supplies the work,
    /// independently of large endpoint storage. It is supplied to/returned from
    /// the physical parameter owner, never deposited as resistor heat. The native
    /// descriptor refuses incompatible ideal clamps and out-of-support charge.
    /// Refusal or cancellation leaves the accepted circuit unchanged.
    pub fn rebind_dielectric_capacitors_candidate(
        &self,
        mut dielectrics: Vec<MaterialDielectricCapacitor>,
        inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialStorageCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        canonical_sources(&mut dielectrics)?;
        if dielectrics.len() != self.dielectric_capacitors.len()
            || dielectrics
                .iter()
                .zip(self.dielectric_capacitors.iter())
                .any(|(after, before)| {
                    after.branch != before.branch
                        || after.area_m2 != before.area_m2
                        || after.gap_m != before.gap_m
                })
        {
            return Err(invalid(
                "dielectric source refresh must retain every original branch and fixed geometry",
            ));
        }
        let mut curves = self.capacitor_curves.to_vec();
        for source in &dielectrics {
            poll(&mut cancelled)?;
            let curve = curves
                .iter_mut()
                .find(|curve| curve.branch == source.branch)
                .ok_or_else(|| invalid("retained dielectric source lost its native charge law"))?;
            *curve = source.binding();
        }
        let mut candidate =
            self.rebind_capacitor_curves_with_work(curves, inputs, &mut cancelled)?;
        candidate.circuit.dielectric_capacitors = dielectrics.into();
        poll(&mut cancelled)?;
        Ok(candidate)
    }
}
