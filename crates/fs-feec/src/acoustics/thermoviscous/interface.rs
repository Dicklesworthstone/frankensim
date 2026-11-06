//! Full vector-shear and thermal continuity between resolved fluid regions.
//!
//! With original weak loads B, each region satisfies A x = f + B lambda.
//! The common P2 traces enforce sum(B^T x)=0, including BOTH tangential
//! velocities and theta=T'/T0. The interface Schur operator is assembled from
//! the retained original region LUs; no pressure-only or inviscid connection
//! replaces the original mixed equations. Total traction balances while each
//! region retains its own P1 pressure space. The same physical trace binding
//! also serves the causal time-domain interface.
use super::*;
use body_force::ThermoviscousBodyForceSource;
use heat_source::ThermoviscousHeatSource;

mod schur;
pub(super) mod trace;
use schur::{ConformingSchur, dual};

fn invalid(message: &'static str) -> ThermoviscousError {
    ThermoviscousError::InvalidInput(message)
}

/// Complete original hierarchical P2 trace in a common geometric ordering.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConformingThermoviscousInterface {
    /// Global P2 DOFs: vertices followed by edges matched by endpoint geometry.
    pub trace_dofs: [Vec<usize>; 2],
    /// Active vector-traction basis, indexed into trace_dofs. Coordinates fixed
    /// on BOTH sides are omitted; a one-sided zero still constrains the other.
    /// Tractions are Galerkin representatives modulo work-null loads, so an
    /// omitted coefficient does not assert zero physical stress.
    pub traction_indices: Vec<usize>,
    /// Analogous active basis for inward conductive flux, dual to theta=T'/T0.
    pub heat_flux_indices: Vec<usize>,
}

/// Additional cold interface work, excluding the two already prepared regions.
#[derive(Clone, Copy, Debug)]
pub struct ThermoviscousInterfaceLimits {
    pub max_trace_dofs: usize,
    /// Physical load columns, complete surface masses and complex Schur arrays.
    pub max_dense_bytes: usize,
    /// Numerical continuity/work tolerance, not a spatial accuracy certificate.
    pub relative_residual_tolerance: f64,
}
impl Default for ThermoviscousInterfaceLimits {
    fn default() -> Self {
        Self {
            max_trace_dofs: 256,
            max_dense_bytes: 64 * 1024 * 1024,
            relative_residual_tolerance: 1e-9,
        }
    }
}

/// Two complete fields, with signed internal transfers accounted once per side.
#[derive(Clone, Debug)]
pub struct ConformingThermoviscousResponse {
    /// Original fields and independent volume-loss quadrature. Each input_w
    /// includes its signed mechanical and conductive interface input below.
    /// Original pressure-drive entries on the internal cut are zero: the cut
    /// is driven by the vector traction and heat flux, not those pressure modes.
    pub regions: [ThermoviscousField; 2],
    /// Original P2 traction on the LEFT [Pa]; the right receives its negative.
    pub traction_pa: Vec<[C64; 3]>,
    /// Original P2 conductive flux into the LEFT [W/m2].
    pub inward_heat_flux_w_m2: Vec<C64>,
    pub mechanical_interface_work_w: [f64; 2],
    /// Acoustic available-energy transfer Re(conj(theta)*q)/2 [W].
    /// This is not the first-order caloric heat-flow phasor.
    pub conductive_interface_work_w: [f64; 2],
    /// All external traction, structural, body-force and heat-source work [W].
    pub external_input_w: f64,
}
impl ConformingThermoviscousResponse {
    pub fn balance_residual_w(&self) -> f64 {
        self.external_input_w
            - self
                .regions
                .iter()
                .map(|r| r.power.dissipated_w())
                .sum::<f64>()
    }
}

/// Original region factors plus one bounded complete-interface Schur factor.
pub struct PreparedConformingThermoviscous<'p, 'a, 'c> {
    regions: [&'p PreparedThermoviscous<'a, 'c>; 2],
    binding: trace::TraceBinding,
    external_pressure: Vec<(usize, usize)>,
    solver: ConformingSchur,
    tolerance: f64,
}

impl<'a, 'c> PreparedThermoviscous<'a, 'c> {
    /// Join matching original face triangulations at this same frequency.
    ///
    /// Complete vector velocity and temperature are continuous. Pressure modes
    /// touching the cut are removed from the external drives and must lie
    /// entirely on it. Molecular complex Cp, volume forces, heat sources and
    /// disjoint moving/conductive walls retain their original owners.
    ///
    /// Fixed geometry and a common equilibrium pressure/temperature are required.
    /// Nonmatching meshes, structural traces overlapping the cut and singular
    /// interface loads refuse. No artificial interface friction is introduced.
    pub fn connect_conforming<'p>(
        &'p self,
        right: &'p Self,
        left_faces: &[[u32; 3]],
        right_faces: &[[u32; 3]],
        limits: ThermoviscousInterfaceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedConformingThermoviscous<'p, 'a, 'c>, ThermoviscousError> {
        poll(cancelled)?;
        if self.frequency_hz != right.frequency_hz
            || self.model.medium.pressure_pa != right.model.medium.pressure_pa
            || self.model.medium.temperature_k != right.model.medium.temperature_k
            || !limits.relative_residual_tolerance.is_finite()
            || limits.relative_residual_tolerance <= 0.0
            || limits.relative_residual_tolerance > 1e-6
        {
            return Err(invalid(
                "common harmonic frequency/equilibrium and bounded positive tolerance required",
            ));
        }
        let regions = [self, right];
        let binding = trace::prepare(
            regions.map(|r| r.model),
            [left_faces, right_faces],
            limits.max_trace_dofs.min(MAX_LNS_DOFS),
            limits.max_dense_bytes,
            cancelled,
        )?;
        let count =
            3 * binding.layout.traction_indices.len() + binding.layout.heat_flux_indices.len();
        let n = self.model.n + right.model.n;
        let bytes = count
            .checked_mul(count)
            .and_then(|v| v.checked_mul(48))
            .and_then(|v| {
                n.checked_mul(count)
                    .and_then(|k| k.checked_mul(8))
                    .and_then(|k| v.checked_add(k))
            })
            .ok_or(ThermoviscousError::NonFinite)?;
        if bytes > limits.max_dense_bytes {
            return Err(ThermoviscousError::Budget {
                what: "conforming harmonic interface bytes",
                requested: bytes,
                limit: limits.max_dense_bytes,
            });
        }
        let solver = ConformingSchur::prepare(
            [&binding.loads[0], &binding.loads[1]],
            regions.map(|source| source.model.n),
            count,
            limits.relative_residual_tolerance,
            &mut |side, rhs, cancelled| solve_region(regions[side], rhs, cancelled),
            cancelled,
        )?;
        let mut external_pressure = Vec::new();
        for side in 0..2 {
            for mode in 0..regions[side].model.loads.len() {
                if !binding.selected[side].contains(&mode) {
                    external_pressure.push((side, mode));
                }
            }
        }
        Ok(PreparedConformingThermoviscous {
            regions,
            binding,
            external_pressure,
            solver,
            tolerance: limits.relative_residual_tolerance,
        })
    }
}

impl PreparedConformingThermoviscous<'_, '_, '_> {
    pub fn interface(&self) -> &ConformingThermoviscousInterface {
        &self.binding.layout
    }
    pub fn external_pressure_modes(&self) -> &[(usize, usize)] {
        &self.external_pressure
    }

    /// Original pressure drives in external_pressure_modes order, and complete
    /// structural-force coordinates per region. No prepared owner is mutated.
    pub fn response(
        &self,
        external_traction_pa: &[C64],
        wall_forces: [&[C64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConformingThermoviscousResponse, ThermoviscousError> {
        self.response_with_volume_sources(
            external_traction_pa,
            wall_forces,
            [&[]; 2],
            [&[]; 2],
            cancelled,
        )
    }

    /// Simultaneous independently phased original sources. Preparation is reused;
    /// source sampling, constitutive evaluation and dense factorization stay cold.
    #[allow(clippy::too_many_arguments)]
    pub fn response_with_volume_sources(
        &self,
        external_traction_pa: &[C64],
        wall_forces: [&[C64]; 2],
        force_sources: [&[(&ThermoviscousBodyForceSource<'_, '_>, C64)]; 2],
        heat_sources: [&[(&ThermoviscousHeatSource<'_, '_>, C64)]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConformingThermoviscousResponse, ThermoviscousError> {
        poll(cancelled)?;
        if external_traction_pa.len() != self.external_pressure.len()
            || external_traction_pa.iter().any(|&z| !finite(z))
        {
            return Err(invalid(
                "complete finite external interface drives required",
            ));
        }
        let count = 3 * self.binding.layout.traction_indices.len()
            + self.binding.layout.heat_flux_indices.len();
        let mechanical = 3 * self.binding.layout.traction_indices.len();
        let mut drive: [Vec<C64>; 2] =
            std::array::from_fn(|s| vec![C64::ZERO; self.regions[s].model.loads.len()]);
        for (&(side, mode), &value) in self.external_pressure.iter().zip(external_traction_pa) {
            drive[side][mode] = value;
        }
        let mut rhs: [Vec<C64>; 2] =
            std::array::from_fn(|s| vec![C64::ZERO; self.regions[s].model.n]);
        for side in 0..2 {
            let model = self.regions[side].model;
            for (load, &pressure) in model.loads.iter().zip(&drive[side]) {
                for (entry, &coefficient) in rhs[side].iter_mut().zip(load) {
                    *entry = *entry + pressure.scale(coefficient);
                }
            }
            heat_source::add_sources(
                model,
                &mut rhs[side],
                wall_forces[side],
                force_sources[side],
                heat_sources[side],
                cancelled,
            )?;
        }
        let (fields, reaction) = self.solver.response(
            [&self.binding.loads[0], &self.binding.loads[1]],
            mechanical,
            self.regions.map(|source| 3 * source.model.velocity_count),
            self.regions.map(|source| source.model.temperature_start),
            [&rhs[0], &rhs[1]],
            &mut |side, rhs, cancelled| solve_region(self.regions[side], rhs, cancelled),
            cancelled,
        )?;
        let mut region_fields = Vec::with_capacity(2);
        let mut mechanical_work = [0.0; 2];
        let mut conductive_work = [0.0; 2];
        let mut external_input_w = 0.0;
        for side in 0..2 {
            let mut field = heat_source::field_with_source_work(
                self.regions[side],
                &drive[side],
                wall_forces[side],
                force_sources[side],
                &fields[side],
                &rhs[side],
                cancelled,
            )?;
            external_input_w += field.power.input_w;
            let mut work = [Dd::ZERO; 2];
            for a in 0..count {
                let flow = dual(&self.binding.loads[side], count, a, &fields[side]);
                let channel = usize::from(a >= mechanical);
                work[channel] = work[channel]
                    + Dd::from_f64(reaction[a].re) * Dd::from_f64(flow.re)
                    + Dd::from_f64(reaction[a].im) * Dd::from_f64(flow.im);
            }
            mechanical_work[side] = 0.5 * work[0].to_f64();
            conductive_work[side] = 0.5 * work[1].to_f64();
            field.power.input_w += mechanical_work[side] + conductive_work[side];
            region_fields.push(field);
        }
        let regions: [ThermoviscousField; 2] =
            region_fields.try_into().expect("two original regions");
        let nv = self.binding.layout.traction_indices.len();
        let response = ConformingThermoviscousResponse {
            regions,
            traction_pa: (0..nv)
                .map(|i| std::array::from_fn(|c| reaction[c * nv + i]))
                .collect(),
            inward_heat_flux_w_m2: reaction[mechanical..].to_vec(),
            mechanical_interface_work_w: mechanical_work,
            conductive_interface_work_w: conductive_work,
            external_input_w,
        };
        let scale = external_input_w.abs()
            + response
                .regions
                .iter()
                .map(|r| r.power.dissipated_w().abs())
                .sum::<f64>()
            + mechanical_work
                .iter()
                .chain(&conductive_work)
                .map(|v| v.abs())
                .sum::<f64>();
        if !scale.is_finite()
            || !response.balance_residual_w().is_finite()
            || response.balance_residual_w().abs() > self.tolerance * scale
            || mechanical_work.iter().sum::<f64>().abs() > self.tolerance * scale
            || conductive_work.iter().sum::<f64>().abs() > self.tolerance * scale
        {
            return Err(invalid(
                "conforming harmonic interface power residual exceeded",
            ));
        }
        poll(cancelled)?;
        Ok(response)
    }
}

fn solve_region(
    source: &PreparedThermoviscous<'_, '_>,
    rhs: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, ThermoviscousError> {
    source.model.solve_system(
        &source.factor,
        &source.rows,
        &source.columns,
        std::f64::consts::TAU * source.frequency_hz,
        source.cp,
        &source.wall_admittances,
        None,
        rhs,
        cancelled,
    )
}

#[cfg(test)]
mod tests;
