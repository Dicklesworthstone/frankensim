//! Accepted reservoir histories and independent transfer-error estimates.

use super::super::super::{
    SpeciesBoundaryDiffusionReport, SpeciesBoundaryTransfer, SpeciesDiffusionBoundary,
};
use super::super::boundary::{
    BoundaryEnthalpyAccount, BoundaryEnthalpySums, SpeciesBoundaryEnthalpyDiffusionReport,
    SpeciesBoundaryEnthalpyTransfer,
};
use super::{
    AdaptiveSpeciesEnthalpyDiffusionControls, AdaptiveSpeciesEnthalpyDiffusionStep, CancelGate,
    ChtError, EnergySum, VoxelDomain, nonnegative, poll_every, positive, refuse,
};

/// Coupled field accuracy plus independent per-face reservoir throughput accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls<const N: usize> {
    pub integrator: AdaptiveSpeciesEnthalpyDiffusionControls<N>,
    /// Each constituent retains its own absolute integrated-transfer allowance.
    pub absolute_boundary_transfer_tolerance_kg: [f64; N],
    /// Sum of per-face signed/positive/negative energy-history estimates [J].
    /// Opposing reservoirs cannot hide behind a stationary endpoint.
    pub absolute_boundary_enthalpy_tolerance_j: f64,
}

/// Original accepted fine trajectory and all external transfers already applied.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionStep<const N: usize, P = ()> {
    pub step: AdaptiveSpeciesEnthalpyDiffusionStep<N, P>,
    /// Integrated transfers cover the complete interval. Face conductances are
    /// the LAST accepted midpoint coefficients, not an interval average.
    pub boundary: SpeciesBoundaryDiffusionReport<N>,
    pub thermal_boundary: SpeciesBoundaryEnthalpyDiffusionReport,
    pub estimated_boundary_transfer_error_kg: [f64; N],
    pub estimated_boundary_enthalpy_error_j: f64,
}

#[derive(Clone, Copy)]
pub(super) struct BoundaryOptions<'a, const N: usize> {
    pub boundaries: &'a [SpeciesDiffusionBoundary<N>],
    pub mass_tolerance: [f64; N],
    pub heat_tolerance: f64,
}

impl<const N: usize> BoundaryOptions<'_, N> {
    pub fn validate(&self) -> Result<(), ChtError> {
        for value in self.mass_tolerance {
            positive(value)?;
        }
        positive(self.heat_tolerance)?;
        Ok(())
    }

    pub fn competing_maps(
        &self,
        domain: &VoxelDomain,
        diffusivity: &[f64],
        gate: &CancelGate,
    ) -> Result<bool, ChtError> {
        super::super::boundary::boundary_maps_compete(
            domain,
            diffusivity,
            self.boundaries,
            gate,
        )
    }

    /// Coarse/current histories, live native face readbacks and final outputs.
    /// Native solve scratch is charged separately by its original storage gate.
    pub fn retained_bytes(&self) -> Option<usize> {
        let per_face = 3usize
            .checked_mul(core::mem::size_of::<HistoryFace<N>>())?
            .checked_add(3usize.checked_mul(core::mem::size_of::<SpeciesBoundaryTransfer<N>>())?)?
            .checked_add(
                3usize.checked_mul(core::mem::size_of::<SpeciesBoundaryEnthalpyTransfer>())?,
            )?;
        self.boundaries
            .len()
            .checked_mul(per_face)?
            .checked_add(3usize.checked_mul(core::mem::size_of::<BoundaryHistory<N>>())?)?
            .checked_add(core::mem::size_of::<
                AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionStep<N>,
            >())?
            .checked_add(2048)
    }
}

struct HistoryFace<const N: usize> {
    cell: usize,
    face: usize,
    conductance: f64,
    sums: BoundaryEnthalpySums<N>,
}

pub(super) struct BoundaryHistory<const N: usize> {
    faces: Vec<HistoryFace<N>>,
    external: EnergySum,
    positive: EnergySum,
    negative: EnergySum,
    mass: [EnergySum; N],
    carrier_defect: EnergySum,
}

fn add(sum: &mut EnergySum, value: EnergySum) -> Result<(), ChtError> {
    sum.add(value.sum);
    sum.add(value.correction);
    if !sum.sum.is_finite() || !sum.correction.is_finite() || !sum.value().is_finite() {
        return Err(refuse("adaptive reservoir history is unrepresentable"));
    }
    Ok(())
}

fn difference(a: EnergySum, b: EnergySum) -> Result<f64, ChtError> {
    let mut result = a;
    result.add(-b.sum);
    result.add(-b.correction);
    let value = result.value().abs();
    nonnegative(value)?;
    Ok(value)
}

fn history_error(a: [EnergySum; 3], b: [EnergySum; 3]) -> Result<f64, ChtError> {
    Ok(difference(a[0], b[0])?
        .max(difference(a[1], b[1])?)
        .max(difference(a[2], b[2])?)
        / 3.0)
}

impl<const N: usize> BoundaryHistory<N> {
    pub fn new(boundaries: &[SpeciesDiffusionBoundary<N>]) -> Self {
        Self {
            faces: boundaries
                .iter()
                .map(|b| HistoryFace {
                    cell: b.cell,
                    face: b.face,
                    conductance: 0.0,
                    sums: BoundaryEnthalpySums::default(),
                })
                .collect(),
            external: EnergySum::default(),
            positive: EnergySum::default(),
            negative: EnergySum::default(),
            mass: [EnergySum::default(); N],
            carrier_defect: EnergySum::default(),
        }
    }

    pub fn accumulate(
        &mut self,
        mass: &SpeciesBoundaryDiffusionReport<N>,
        account: &BoundaryEnthalpyAccount<N>,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        if self.faces.len() != account.faces.len() || self.faces.len() != mass.faces.len() {
            return Err(refuse(
                "adaptive reservoir history lost its canonical faces",
            ));
        }
        for (i, ((out, source), original)) in self
            .faces
            .iter_mut()
            .zip(&account.faces)
            .zip(&mass.faces)
            .enumerate()
        {
            poll_every(gate, i)?;
            if (out.cell, out.face) != (original.cell, original.face) {
                return Err(refuse("adaptive reservoir face identity changed"));
            }
            out.conductance = original.conductance_kg_s;
            add(&mut out.sums.signed, source.signed)?;
            add(&mut out.sums.positive, source.positive)?;
            add(&mut out.sums.negative, source.negative)?;
            for k in 0..N {
                add(&mut out.sums.mass_signed[k], source.mass_signed[k])?;
                add(&mut out.sums.mass_incoming[k], source.mass_incoming[k])?;
                add(&mut out.sums.mass_outgoing[k], source.mass_outgoing[k])?;
            }
        }
        add(&mut self.external, account.external)?;
        add(&mut self.positive, account.positive)?;
        add(&mut self.negative, account.negative)?;
        for k in 0..N {
            add(&mut self.mass[k], account.mass_into_domain[k])?;
        }
        self.carrier_defect.add(mass.absolute_carrier_transfer_kg);
        nonnegative(self.carrier_defect.value())?;
        Ok(())
    }

    pub fn external_energy(&self) -> EnergySum {
        self.external
    }
    pub fn external_mass(&self) -> [EnergySum; N] {
        self.mass
    }

    pub fn compare(
        &self,
        fine: &Self,
        controls: BoundaryOptions<'_, N>,
        gate: &CancelGate,
    ) -> Result<([f64; N], f64, f64), ChtError> {
        if self.faces.len() != fine.faces.len() {
            return Err(refuse("adaptive reservoir comparison lost a face"));
        }
        let mut mass = [EnergySum::default(); N];
        let mut heat = EnergySum::default();
        for (i, (a, b)) in self.faces.iter().zip(&fine.faces).enumerate() {
            poll_every(gate, i)?;
            if (a.cell, a.face) != (b.cell, b.face) {
                return Err(refuse(
                    "adaptive reservoir comparison changed face identity",
                ));
            }
            heat.add(history_error(
                [a.sums.signed, a.sums.positive, a.sums.negative],
                [b.sums.signed, b.sums.positive, b.sums.negative],
            )?);
            for (k, sum) in mass.iter_mut().enumerate() {
                sum.add(history_error(
                    [
                        a.sums.mass_signed[k],
                        a.sums.mass_incoming[k],
                        a.sums.mass_outgoing[k],
                    ],
                    [
                        b.sums.mass_signed[k],
                        b.sums.mass_incoming[k],
                        b.sums.mass_outgoing[k],
                    ],
                )?);
            }
        }
        let mass = mass.map(EnergySum::value);
        let heat = heat.value();
        nonnegative(heat)?;
        let mut ratio = heat / controls.heat_tolerance;
        for k in 0..N {
            nonnegative(mass[k])?;
            ratio = ratio.max(mass[k] / controls.mass_tolerance[k]);
        }
        // Finite error divided by a very small positive tolerance may be
        // infinite. It remains an unresolved refinement ratio, not a bad law.
        Ok((mass, heat, ratio))
    }

    pub fn finish(
        self,
        gate: &CancelGate,
    ) -> Result<
        (
            SpeciesBoundaryDiffusionReport<N>,
            SpeciesBoundaryEnthalpyDiffusionReport,
        ),
        ChtError,
    > {
        let mut mass_faces = Vec::with_capacity(self.faces.len());
        let mut heat_faces = Vec::with_capacity(self.faces.len());
        for (i, entry) in self.faces.into_iter().enumerate() {
            poll_every(gate, i)?;
            mass_faces.push(SpeciesBoundaryTransfer {
                cell: entry.cell,
                face: entry.face,
                conductance_kg_s: entry.conductance,
                mass_into_cell_kg: entry.sums.mass_signed.map(EnergySum::value),
                incoming_mass_kg: entry.sums.mass_incoming.map(EnergySum::value),
                outgoing_mass_kg: entry.sums.mass_outgoing.map(EnergySum::value),
            });
            heat_faces.push(SpeciesBoundaryEnthalpyTransfer {
                cell: entry.cell,
                face: entry.face,
                enthalpy_into_cell_j: entry.sums.signed.value(),
                positive_enthalpy_into_cell_j: entry.sums.positive.value(),
                negative_enthalpy_into_cell_j: entry.sums.negative.value(),
            });
        }
        Ok((
            SpeciesBoundaryDiffusionReport {
                faces: mass_faces,
                mass_into_domain_kg: self.mass.map(EnergySum::value),
                absolute_carrier_transfer_kg: self.carrier_defect.value(),
            },
            SpeciesBoundaryEnthalpyDiffusionReport {
                faces: heat_faces,
                enthalpy_into_domain_j: self.external.value(),
                positive_enthalpy_into_domain_j: self.positive.value(),
                negative_enthalpy_into_domain_j: self.negative.value(),
            },
        ))
    }
}
