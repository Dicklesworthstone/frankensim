//! Full-stick vector tractions from the original mixed-mode JKR equilibrium.
//!
//! Mindlin (1949), Compliance of Elastic Bodies in Contact, equation 76:
//! q_i(r) = T_i / (2 pi a sqrt(a*a-r*r)). This is the SAME circular,
//! elastically decoupled, full-stick model used by SourceJkrMixedMode. Its
//! singular shear field has the tensile flat-punch radial centroid, so force,
//! first moments and work are integrated exactly without rim sampling.
//!
//! Ordinary JKR/Maugis profiles have zero shear. The original normal-only APIs
//! remain unchanged; these vector operations are an explicit additional choice.
//! Partial slip, torsional compliance and finite-body contact equilibrium are
//! separate physical owners. A prescribed traction does not establish them.

use super::super::{AdhesiveTractionError, AdhesiveTractionProfile};
use super::{
    MAX_POLAR_TRACTION_SECTORS, PolarRingScratch, PolarTractionTransfer, PolarTransferBudget,
    RadialMomentBudget, SectorAdhesiveTraction, SectorTractionError, angles, check_budget,
    check_sector, edge_root, finite, poll, sector_geometry, sector_with_geometry,
};
use core::f64::consts::TAU;
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// A complete local traction sector. Local component order is x shear, y
/// shear, normal; moments use position relative to the original disk centre.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SectorAdhesiveVectorTraction {
    /// Original normal pressure integral and its quadrature estimate.
    pub normal: SectorAdhesiveTraction,
    pub tangential_force_n: [f64; 2],
    /// Entry [i][j] is integral of position_j * shear_traction_i dA, N m.
    /// Rows are x/y shear components; columns are x/y position coordinates.
    pub tangential_first_moment_n_m: [[f64; 2]; 2],
}

impl SectorAdhesiveVectorTraction {
    /// Integrated local force [T_x,T_y,N], N.
    #[must_use]
    pub fn force_n(&self) -> [f64; 3] {
        [
            self.tangential_force_n[0],
            self.tangential_force_n[1],
            self.normal.normal_force_n,
        ]
    }

    /// Entry [i][j] is integral of position_j * traction_i dA, N m.
    #[must_use]
    pub fn first_moments_n_m(&self) -> [[f64; 2]; 3] {
        [
            self.tangential_first_moment_n_m[0],
            self.tangential_first_moment_n_m[1],
            self.normal.first_moment_n_m,
        ]
    }

    /// Force and torque about the original centre in the right-handed local
    /// frame. Opposite bodies receive opposite wrenches.
    #[must_use]
    pub fn wrench(&self) -> ([f64; 3], [f64; 3]) {
        let m = self.first_moments_n_m();
        (
            self.force_n(),
            [
                m[2][1],
                -m[2][0],
                (Dd::from_f64(m[1][0]) - Dd::from_f64(m[0][1])).to_f64(),
            ],
        )
    }

    /// Work-conjugate power against v_i(x,y)=translation_i+gradient_ij*x_j.
    /// Components are [x shear,y shear,normal]; gradient columns are [x,y].
    /// The error estimate covers only the original normal-moment quadrature.
    /// Full-stick shear integrals are analytic; rounding/model error is excluded.
    pub fn affine_power(
        &self,
        translation_m_per_s: [f64; 3],
        gradient_per_s: [[f64; 2]; 3],
    ) -> Result<(f64, f64), SectorTractionError> {
        let force = self.force_n();
        let moment = self.first_moments_n_m();
        let mut power = Dd::ZERO;
        for i in 0..3 {
            finite(translation_m_per_s[i])?;
            power = power + Dd::from_f64(force[i]) * Dd::from_f64(translation_m_per_s[i]);
            for j in 0..2 {
                finite(gradient_per_s[i][j])?;
                power = power + Dd::from_f64(moment[i][j]) * Dd::from_f64(gradient_per_s[i][j]);
            }
        }
        Ok((
            finite(power.to_f64())?,
            finite(
                self.normal.moment_error_estimate_n_m
                    * gradient_per_s[2][0].hypot(gradient_per_s[2][1]),
            )?,
        ))
    }
}

/// Caller-owned bounded scratch; a refusal may change scratch but never output.
#[derive(Clone, Copy, Debug, Default)]
pub struct VectorPolarRingScratch {
    normal: PolarRingScratch,
    tangential_force_n: [f64; 2],
    tangential_radial_moment_n_m: [f64; 2],
}

/// Original source identities and all signed resultant residuals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VectorPolarTractionTransfer {
    pub normal: PolarTractionTransfer,
    pub tangential_force_residual_n: [f64; 2],
    /// Rows are x/y shear components; columns are x/y position coordinates.
    pub tangential_first_moment_residual_n_m: [[f64; 2]; 2],
}

impl AdhesiveTractionProfile {
    /// Pointwise local traction [T_x,T_y,N], Pa. At a JKR edge use integrated
    /// sectors instead: both ideal tensile and full-stick shear fields diverge.
    pub fn traction_pa(&self, radius_m: f64) -> Result<[f64; 3], AdhesiveTractionError> {
        let normal = self.pressure_pa(radius_m)?;
        let mut tangential = [0.0; 2];
        if radius_m < self.contact_radius_m() && self.tangential_force_n().iter().any(|&q| q != 0.0)
        {
            let a = self.contact_radius_m();
            let root = edge_root(a, radius_m);
            for (out, force) in tangential.iter_mut().zip(self.tangential_force_n()) {
                *out = finite(((force / a) / TAU) / (a * root))?;
            }
        }
        Ok([tangential[0], tangential[1], normal])
    }

    /// Integrate one sector of the same normal field and full-stick shear
    /// field. Exact angular factors and cached radial centroids carry couples
    /// even when a sector's signed normal force is zero.
    #[allow(clippy::too_many_arguments)]
    pub fn vector_sector(
        &self,
        inner_m: f64,
        outer_m: f64,
        start_rad: f64,
        end_rad: f64,
        budget: RadialMomentBudget,
        gate: &CancelGate,
    ) -> Result<SectorAdhesiveVectorTraction, SectorTractionError> {
        check_budget(budget)?;
        angles(start_rad, end_rad)?;
        poll(gate)?;
        let mut evaluations = 0;
        let radial = self.vector_radial_load(
            inner_m,
            outer_m,
            budget,
            &mut evaluations,
            gate,
            end_rad - start_rad != TAU,
        )?;
        let cell = vector_sector_load(radial, start_rad, end_rad);
        check_vector_sector(cell)?;
        poll(gate)?;
        Ok(cell)
    }

    /// Transfer a complete radial-major/angular-minor grid. Normal and shear
    /// forces and every first moment pass the original absolute vector-norm
    /// gates before publication. One radial calculation is shared by all angles;
    /// the full-stick field adds no quadrature and no per-cell allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn transfer_vector_polar_grid(
        &self,
        radii_m: &[f64],
        angles_rad: &[f64],
        output: &mut [SectorAdhesiveVectorTraction],
        scratch: &mut [VectorPolarRingScratch],
        budget: PolarTransferBudget,
        gate: &CancelGate,
    ) -> Result<VectorPolarTractionTransfer, SectorTractionError> {
        check_budget(budget.radial_moment)?;
        poll(gate)?;
        if radii_m.len() < 2
            || angles_rad.len() < 2
            || radii_m[0] != 0.0
            || angles_rad[0] != 0.0
            || *angles_rad.last().unwrap() != TAU
            || *radii_m.last().unwrap() < self.interaction_radius_m()
            || scratch.len() != radii_m.len() - 1
            || (radii_m.len() - 1).checked_mul(angles_rad.len() - 1) != Some(output.len())
            || output.len() > MAX_POLAR_TRACTION_SECTORS
            || !budget.force_absolute_tolerance_n.is_finite()
            || budget.force_absolute_tolerance_n < 0.0
            || !budget.moment_absolute_tolerance_n_m.is_finite()
            || budget.moment_absolute_tolerance_n_m < 0.0
        {
            return Err(SectorTractionError::Invalid(
                "complete vector polar grid, storage and resultant gates",
            ));
        }
        for edge in angles_rad.windows(2) {
            angles(edge[0], edge[1])?;
        }
        let mut evaluations = 0;
        for (edge, staged) in radii_m.windows(2).zip(scratch.iter_mut()) {
            poll(gate)?;
            *staged = self.vector_radial_load(
                edge[0],
                edge[1],
                budget.radial_moment,
                &mut evaluations,
                gate,
                angles_rad.len() > 2,
            )?;
        }
        let mut forces = [Dd::ZERO; 3];
        let mut moments = [[Dd::ZERO; 2]; 3];
        let mut largest_error = 0.0_f64;
        for radial in scratch.iter().copied() {
            poll(gate)?;
            for edge in angles_rad.windows(2) {
                let cell = vector_sector_load(radial, edge[0], edge[1]);
                check_vector_sector(cell)?;
                let force = cell.force_n();
                let moment = cell.first_moments_n_m();
                for i in 0..3 {
                    forces[i] = forces[i] + Dd::from_f64(force[i]);
                    for j in 0..2 {
                        moments[i][j] = moments[i][j] + Dd::from_f64(moment[i][j]);
                    }
                }
                largest_error = largest_error.max(cell.normal.moment_error_estimate_n_m);
            }
        }
        let original = [
            self.tangential_force_n()[0],
            self.tangential_force_n()[1],
            self.normal_force_n(),
        ];
        let mut force_residual_n = [0.0; 3];
        let mut first_moment_residual_n_m = [[0.0; 2]; 3];
        let mut force_norm = 0.0_f64;
        let mut moment_norm = 0.0_f64;
        for i in 0..3 {
            force_residual_n[i] = finite((forces[i] - Dd::from_f64(original[i])).to_f64())?;
            force_norm = force_norm.hypot(force_residual_n[i]);
            for j in 0..2 {
                first_moment_residual_n_m[i][j] = finite(moments[i][j].to_f64())?;
                moment_norm = moment_norm.hypot(first_moment_residual_n_m[i][j]);
            }
        }
        if force_norm > budget.force_absolute_tolerance_n
            || moment_norm > budget.moment_absolute_tolerance_n_m
        {
            return Err(SectorTractionError::VectorResultant {
                force_residual_n,
                first_moment_residual_n_m,
            });
        }
        let normal = PolarTractionTransfer {
            model_identity: self.model_identity(),
            normal_state_identity: self.normal_state_identity(),
            authority: self.authority(),
            force_residual_n: force_residual_n[2],
            first_moment_residual_n_m: first_moment_residual_n_m[2],
            maximum_sector_moment_error_estimate_n_m: largest_error,
            radial_evaluations: evaluations,
        };
        poll(gate)?;
        let sectors = angles_rad.len() - 1;
        for (radial, cells) in scratch.iter().copied().zip(output.chunks_mut(sectors)) {
            for (edge, cell) in angles_rad.windows(2).zip(cells) {
                *cell = vector_sector_load(radial, edge[0], edge[1]);
            }
        }
        Ok(VectorPolarTractionTransfer {
            normal,
            tangential_force_residual_n: [force_residual_n[0], force_residual_n[1]],
            tangential_first_moment_residual_n_m: [
                first_moment_residual_n_m[0],
                first_moment_residual_n_m[1],
            ],
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn vector_radial_load(
        &self,
        inner: f64,
        outer: f64,
        budget: RadialMomentBudget,
        evaluations: &mut usize,
        gate: &CancelGate,
        need_moments: bool,
    ) -> Result<VectorPolarRingScratch, SectorTractionError> {
        let normal = if need_moments {
            self.radial_load(inner, outer, budget, evaluations, gate)?
        } else {
            PolarRingScratch {
                ring: self.annulus(inner, outer)?,
                ..PolarRingScratch::default()
            }
        };
        let mut radial = VectorPolarRingScratch {
            normal,
            ..VectorPolarRingScratch::default()
        };
        let a = self.contact_radius_m();
        let lo = inner.min(a);
        let hi = outer.min(a);
        if hi > lo && self.tangential_force_n().iter().any(|&q| q != 0.0) {
            // sqrt(1-lo²/a²)-sqrt(1-hi²/a²), factored for thin rings.
            let weight = finite(
                ((hi - lo) / a) * ((hi / a + lo / a) / (edge_root(a, lo) + edge_root(a, hi))),
            )?;
            for i in 0..2 {
                radial.tangential_force_n[i] = finite(self.tangential_force_n()[i] * weight)?;
                radial.tangential_radial_moment_n_m[i] =
                    finite(radial.tangential_force_n[i] * (normal.tensile_centroid_m / TAU))?;
            }
        }
        Ok(radial)
    }
}

fn vector_sector_load(
    radial: VectorPolarRingScratch,
    start: f64,
    end: f64,
) -> SectorAdhesiveVectorTraction {
    let (fraction, angular) = sector_geometry(start, end);
    SectorAdhesiveVectorTraction {
        normal: sector_with_geometry(radial.normal, fraction, angular),
        tangential_force_n: radial.tangential_force_n.map(|q| q * fraction),
        tangential_first_moment_n_m: radial
            .tangential_radial_moment_n_m
            .map(|q| [q * angular[0], q * angular[1]]),
    }
}

fn check_vector_sector(cell: SectorAdhesiveVectorTraction) -> Result<(), SectorTractionError> {
    check_sector(cell.normal)?;
    for value in cell
        .tangential_force_n
        .into_iter()
        .chain(cell.tangential_first_moment_n_m.into_iter().flatten())
    {
        finite(value)?;
    }
    Ok(())
}
