//! Complete matching P2 traces of the original resolved equations.
//!
//! A Cartesian traction t [Pa] is applied to the left and -t to the right.
//! Conductive flux q [W/m2] is inward to the left and outward from the right.
//! The original surface mass makes their duals integral(phi v dA) and
//! integral(phi theta dA), theta=T'/T0. Thus the common-flow constraint is
//! equality of the COMPLETE active P2 velocity and temperature traces. The
//! heat term in the acoustic work ledger is theta*q, not the caloric flux q.
use super::*;

pub use super::super::super::interface::ConformingThermoviscousInterface;

/// Same complete trajectory and original work owners as the normal connection.
pub struct PreparedConformingThermoviscousTimePair<'p, 'a, 'c> {
    inner: PreparedThermoviscousTimePair<'p, 'a, 'c>,
    layout: ConformingThermoviscousInterface,
    heat_flux_areas_m2: Vec<f64>,
}

/// Fully admitted fields and separately integrated physical interface transfers.
#[derive(Clone, Debug)]
pub struct ConformingThermoviscousTimePairStep {
    pub regions: [ThermoviscousTimeStep; 2],
    /// Integral(t.v dA dt), inward to each region [J].
    pub mechanical_interface_work_j: [f64; 2],
    /// Integral((T'/T0) q_in dA dt) [J], the signed transfer of acoustic
    /// available energy. This is distinct from actual caloric heat below.
    pub conductive_interface_work_j: [f64; 2],
    pub external_pressure_work_j: f64,
    /// P2 vector traction ON THE LEFT [Pa], in `traction_indices` order.
    /// The traction on the right is its negative in the same Cartesian frame.
    /// On a constrained trace this is the active-space Galerkin representative
    /// modulo work-null loads; omitted coefficients are not physical zero stress.
    pub mean_traction_pa: Vec<[f64; 3]>,
    pub endpoint_traction_pa: Vec<[f64; 3]>,
    /// P2 conductive flux INWARD TO THE LEFT [W/m2], `heat_flux_indices` order.
    /// On a constrained thermal trace this is only the active-space Galerkin
    /// representative. Unobservable flux through fixed-temperature DOFs is not
    /// inferred or silently assigned zero.
    pub mean_inward_heat_flux_w_m2: Vec<f64>,
    pub endpoint_inward_heat_flux_w_m2: Vec<f64>,
    /// Integral(q_in dA dt) [J] for each region, only when the COMPLETE thermal
    /// trace is active. Otherwise a work-null load can have a nonzero integral
    /// and the actual caloric transfer is unidentifiable (`None`). Identified
    /// transfers cancel and are NOT added to the quadratic acoustic ledger.
    pub caloric_heat_transfer_j: Option<[f64; 2]>,
}
impl ConformingThermoviscousTimePairStep {
    pub fn balance_residual_j(&self) -> f64 {
        self.regions
            .iter()
            .map(|r| {
                r.work.delta_energy_j + r.work.dissipated_j()
                    - r.work.structural_work_j
                    - r.work.body_force_work_j
                    - r.work.heat_source_work_j
            })
            .sum::<f64>()
            - self.external_pressure_work_j
    }
}

impl<'a, 'c> PreparedThermoviscousTime<'a, 'c> {
    /// Join matching original face triangulations with complete P2 vector and
    /// temperature continuity. Every pressure mode touching the cut is removed
    /// from the external inputs; such a mode must lie entirely on the cut.
    /// Other pressure and structural-force coordinates retain their ordering.
    ///
    /// Original stationary constraints at the edge are retained, including a
    /// constraint present on only one side. A moving structural trace sharing
    /// a cut DOF is refused: its kinematic constraint must have one owner.
    /// Pressure remains in each region's original P1 volume space; it is not a
    /// new boundary constraint. Total traction, rather than pressure alone,
    /// balances across a viscous material interface. No mean-flow, geometry,
    /// species transport, mean thermal reservoir or contact resistance is added.
    pub fn connect_conforming<'p>(
        &'p self,
        right: &'p Self,
        left_faces: &[[u32; 3]],
        right_faces: &[[u32; 3]],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedConformingThermoviscousTimePair<'p, 'a, 'c>, ThermoviscousTimeError> {
        let regions = [self, right];
        let n = self.system.state_dim() + right.system.state_dim();
        if n > limits.max_states || n >= MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2 {
            return Err(invalid(
                "complete conforming interface state capacity required",
            ));
        }
        let (layout, selected, columns, heat_flux_areas_m2) =
            binding(regions, [left_faces, right_faces], limits, cancelled)?;
        let inner = PreparedThermoviscousTimePair::prepare(
            regions,
            selected,
            Vec::new(),
            Some(columns),
            limits,
            cancelled,
        )?;
        Ok(PreparedConformingThermoviscousTimePair {
            inner,
            layout,
            heat_flux_areas_m2,
        })
    }
}

impl<'p, 'a, 'c> PreparedConformingThermoviscousTimePair<'p, 'a, 'c> {
    pub fn interface(&self) -> &ConformingThermoviscousInterface {
        &self.layout
    }
    pub fn system(&self) -> &PortHamiltonian {
        self.inner.system()
    }
    pub fn duration_s(&self) -> f64 {
        self.inner.duration_s()
    }
    pub fn external_pressure_modes(&self) -> &[(usize, usize)] {
        self.inner.external_pressure_modes()
    }
    pub fn zero_states(&self) -> [ThermoviscousTimeState; 2] {
        self.inner.zero_states()
    }

    /// Both original states remain untouched on every refusal or cancellation.
    /// Initial trace mismatch is rejected; it is never projected away.
    pub fn step(
        &self,
        initial: [&ThermoviscousTimeState; 2],
        external_traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConformingThermoviscousTimePairStep, ThermoviscousTimeError> {
        self.step_with_volume_sources(
            initial,
            external_traction_pa,
            wall_forces_n_per_sqrt_kg,
            [&[]; 2],
            [&[]; 2],
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::acoustics::thermoviscous::time) fn step_with_volume_sources(
        &self,
        initial: [&ThermoviscousTimeState; 2],
        external_traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: [&[f64]; 2],
        force_drive_n: [&[f64]; 2],
        heat_drive_w: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConformingThermoviscousTimePairStep, ThermoviscousTimeError> {
        let frame = self.inner.evaluate_with_volume_sources(
            initial,
            external_traction_pa,
            wall_forces_n_per_sqrt_kg,
            force_drive_n,
            heat_drive_w,
            cancelled,
        )?;
        let nv = self.layout.traction_indices.len();
        let vectors = |reaction: &[f64]| {
            (0..nv)
                .map(|i| std::array::from_fn(|c| reaction[c * nv + i]))
                .collect::<Vec<_>>()
        };
        let mean_traction_pa = vectors(&frame.mean_reaction);
        let endpoint_traction_pa = vectors(&frame.endpoint_reaction);
        let mean_inward_heat_flux_w_m2 = frame.mean_reaction[3 * nv..].to_vec();
        let endpoint_inward_heat_flux_w_m2 = frame.endpoint_reaction[3 * nv..].to_vec();
        let caloric_heat_transfer_j =
            if self.layout.heat_flux_indices.len() == self.layout.trace_dofs[0].len() {
                let transfer: f64 = mean_inward_heat_flux_w_m2
                    .iter()
                    .zip(&self.heat_flux_areas_m2)
                    .map(|(q, area)| q * area * self.duration_s())
                    .sum();
                if !transfer.is_finite() {
                    return Err(invalid(
                        "original conforming caloric heat transfer is nonfinite",
                    ));
                }
                Some([transfer, -transfer])
            } else {
                None
            };
        poll(cancelled)?;
        Ok(ConformingThermoviscousTimePairStep {
            regions: frame.regions,
            mechanical_interface_work_j: frame.mechanical_work_j,
            conductive_interface_work_j: frame.conductive_work_j,
            external_pressure_work_j: frame.external_pressure_work_j,
            mean_traction_pa,
            endpoint_traction_pa,
            mean_inward_heat_flux_w_m2,
            endpoint_inward_heat_flux_w_m2,
            caloric_heat_transfer_j,
        })
    }
}

fn binding(
    regions: [&PreparedThermoviscousTime<'_, '_>; 2],
    faces: [&[[u32; 3]]; 2],
    limits: ThermoviscousTimeLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<
    (
        ConformingThermoviscousInterface,
        [Vec<usize>; 2],
        InterfaceColumns,
        Vec<f64>,
    ),
    ThermoviscousTimeError,
> {
    let original = super::super::super::interface::trace::prepare(
        regions.map(|source| source.model),
        faces,
        limits.max_states,
        limits.preparation.max_dense_bytes,
        cancelled,
    )?;
    let mechanical = 3 * original.layout.traction_indices.len();
    let count = mechanical + original.layout.heat_flux_indices.len();
    let sizes = regions.map(|source| source.system.state_dim());
    let starts = [0, sizes[0]];
    let n = sizes.iter().sum::<usize>();
    let mut loads = vec![0.0; n * count];
    for side in 0..2 {
        let source = regions[side];
        let size = sizes[side];
        let mut physical = vec![0.0; size * count];
        for row in 0..source.model.n {
            poll(cancelled)?;
            for column in 0..count {
                physical[row * count + column] =
                    original.loads[side][row * count + column] / source.scale[row];
            }
        }
        for row in 0..size {
            poll(cancelled)?;
            for column in 0..count {
                loads[(starts[side] + row) * count + column] = (0..size)
                    .map(|k| source.inverse_mass[row * size + k] * physical[k * count + column])
                    .sum();
            }
        }
    }
    Ok((
        original.layout,
        original.selected,
        InterfaceColumns {
            count,
            mechanical,
            loads,
        },
        original.heat_flux_areas_m2,
    ))
}

#[cfg(test)]
mod tests;
