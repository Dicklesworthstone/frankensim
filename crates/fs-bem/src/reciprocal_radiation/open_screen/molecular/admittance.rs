//! A passive physical sheet law in the original relaxing-gas screen.
//!
//! For exp(-i omega t), v_total = v_drive - Y*[p] and the weak equation is
//! (H(k) + i*rho*omega*M_Y)*[p] = i*rho*omega*G^T*rate.
//! Y is measured in m/(Pa*s), not normalized with a complex sound speed.
//! The same P1 interpolation and quadrature supply M_Y and irreversible heat.
//! Bulk molecular relaxation stays entirely in the original material and k.
use super::*;

/// Peak-phasor cycle-average work on the two-sided sheet, in watts.
/// This is a physical channel decomposition, not an independent error bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MolecularScreenPower {
    /// Mechanical work supplied against the fluid reaction. Under illumination
    /// it can be negative because incident sound can drive the structure.
    pub mechanical_power_w: f64,
    /// Irreversible local sheet work, integral Re(Y)*|[p]|^2/2.
    pub sheet_absorbed_power_w: f64,
    /// Signed net work from the sheet into the surrounding fluid:
    /// mechanical_power_w - sheet_absorbed_power_w. This includes molecular
    /// absorption and energy flux, without treating either as an added load.
    /// Incident sources contribute their own power outside this surface.
    /// This balance difference does not bound relative error at cancellation.
    pub fluid_power_w: f64,
}

impl MolecularScreenPower {
    pub(super) fn new(
        mechanical_power_w: f64,
        sheet_absorbed_power_w: f64,
    ) -> Result<Self, HelmholtzError> {
        let fluid_power_w = mechanical_power_w - sheet_absorbed_power_w;
        if !mechanical_power_w.is_finite()
            || !sheet_absorbed_power_w.is_finite()
            || sheet_absorbed_power_w < 0.0
            || !fluid_power_w.is_finite()
        {
            return Err(bad("molecular sheet work is unrepresentable"));
        }
        Ok(Self {
            mechanical_power_w,
            sheet_absorbed_power_w,
            fluid_power_w,
        })
    }
}

impl<'s> PreparedScreen<'s> {
    /// Complete original source work plus one local P1 sheet assembly.
    /// The sheet term is assembled once at cached triangle quadrature points,
    /// with no dense column-by-column application of a diagonal approximation.
    pub fn molecular_velocity_port_work_with_local_admittance(
        &self,
        ports: usize,
    ) -> Result<u64, HelmholtzError> {
        self.molecular_port_boundary_work(ports, true)
    }

    /// Physical velocity ports through a passive porous or locally reacting
    /// sheet, in the SAME relaxing gas on both sides. Every original vertex,
    /// including free-edge vertices, has one Y [m/(Pa*s)] with Re(Y)>=0.
    /// Imaginary Y retains storage; no real/complex sound-speed normalization
    /// is introduced. The caller owns this sheet law's material, frequency
    /// dependence and homogenization validity.
    ///
    /// One H(k)+i*rho*omega*M_Y and LU serve all original source columns.
    /// Normal motion and the pressure jump already represent both faces.
    /// The supplied Y is copied and remains bound to the material/frequency.
    /// An all-zero law preserves the original rigid field arithmetic exactly.
    ///
    /// maximum_jump_entries includes the additional nodal Y row. The original
    /// workspace bound additionally includes that row and one coherent-jump
    /// heat scratch. Active, nonlocal, independently coated face, cavity-backed
    /// and mean-flow laws are outside this local passive sheet model.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_molecular_velocity_ports_with_local_admittance<'p>(
        &'p self,
        source: &MolecularExterior,
        omega: f64,
        rows: &[&[[f64; 3]]],
        admittance_m_per_pa_s: &[C64],
        limits: MolecularScreenPortLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularScreenVelocityPortSolution<'p, 's>, MolecularError> {
        self.solve_molecular_ports_with(
            source,
            omega,
            rows,
            limits,
            Some(admittance_m_per_pa_s),
            None,
            |_, _, _, _, _| Ok(()),
            cancelled,
        )
        .map(|(motion, ())| motion)
    }

    pub(super) fn admit_molecular_local_admittance(
        &self,
        admittance: Option<&[C64]>,
    ) -> Result<(), HelmholtzError> {
        if let Some(values) = admittance {
            if values.len() != self.surface.vertices.len()
                || values.iter().any(|&y| !finite(y) || y.re < 0.0)
            {
                return Err(bad(
                    "molecular sheet needs one finite passive physical admittance per vertex",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn molecular_port_boundary_work(
        &self,
        ports: usize,
        with_sheet: bool,
    ) -> Result<u64, HelmholtzError> {
        let base = u128::from(self.molecular_velocity_port_work(ports)?);
        let extra = if with_sheet {
            let samples = self
                .samples
                .iter()
                .map(|row| row.len() as u128)
                .sum::<u128>();
            128 * samples + 32 * self.surface.vertices.len() as u128
        } else {
            0
        };
        u64::try_from(base + extra).map_err(|_| bad("molecular sheet work overflows"))
    }

    pub(super) fn molecular_port_boundary_payload(
        &self,
        ports: usize,
        with_sheet: bool,
    ) -> Result<u128, HelmholtzError> {
        let extra = if with_sheet {
            32 * self.surface.vertices.len() as u128
        } else {
            0
        };
        Ok(self.molecular_port_payload(ports)? + extra)
    }

    pub(super) fn add_molecular_sheet(
        &self,
        h: &mut [C64],
        rho_omega: f64,
        admittance: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HelmholtzError> {
        if admittance.iter().all(|&y| y == C64::ZERO) {
            return Ok(());
        }
        // Integral Y(x)*phi_i(x)*phi_j(x) is a symmetric cubic P1 form.
        // M*diag(Y) would destroy its work identity for spatially varying Y.
        for (facet, samples) in self.facets.iter().zip(&self.samples) {
            poll(cancelled)?;
            for sample in samples {
                let y = (0..3).fold(C64::ZERO, |sum, a| {
                    sum + admittance[facet.ids[a]].scale(sample.shape[a])
                });
                let weighted = C64::new(-rho_omega * y.im, rho_omega * y.re).scale(sample.weight);
                if !finite(weighted) {
                    return Err(bad("molecular sheet weak coefficient overflows"));
                }
                for a in 0..3 {
                    if let Some(i) = self.dof[facet.ids[a]] {
                        for b in 0..3 {
                            if let Some(j) = self.dof[facet.ids[b]] {
                                let index = i * self.n + j;
                                h[index] =
                                    h[index] + weighted.scale(sample.shape[a] * sample.shape[b]);
                                if !finite(h[index]) {
                                    return Err(bad("molecular sheet weak operator overflows"));
                                }
                            }
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        Ok(())
    }
}

impl MolecularScreenVelocityPortSolution<'_, '_> {
    /// The copied physical nodal Y [m/(Pa*s)] at this source frequency.
    /// None denotes the rigid API; a supplied all-zero row remains present.
    pub fn local_admittance_m_per_pa_s(&self) -> Option<&[C64]> {
        self.local_admittance.as_deref()
    }

    /// Sheet heat of the actual coherent motion field. Original complex rates
    /// combine before squaring; summing unit-column heat would lose interference.
    pub fn sheet_absorbed_power_w(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, HelmholtzError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        if self.local_admittance.is_none() {
            return Ok(0.0);
        }
        let jump = self.pressure_jump(rates, &mut cancelled)?;
        self.sheet_power_from_jump(&jump, &mut cancelled)
    }

    /// Separate mechanical work, local irreversible sheet heat, and signed
    /// net work into the exterior gas. The latter retains molecular loss and
    /// outgoing energy together; it is not a far-field radiation power.
    pub fn power(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularScreenPower, HelmholtzError> {
        let mechanical = self.boundary_work_w(rates, &mut cancelled)?;
        let sheet = self.sheet_absorbed_power_w(rates, &mut cancelled)?;
        MolecularScreenPower::new(mechanical, sheet)
    }

    pub(super) fn sheet_power_from_jump(
        &self,
        jump: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, HelmholtzError> {
        poll(cancelled)?;
        let Some(admittance) = self.local_admittance.as_deref() else {
            return Ok(0.0);
        };
        // Direct conductance work preserves zero loss and very high Q:
        // never form Y*p and cancel its rounded reactive current against p.
        let heat = 0.5 * self.prepared.local_wall_work(admittance, jump, cancelled)?;
        if !heat.is_finite() || heat < 0.0 {
            return Err(bad(
                "coherent molecular sheet absorption is unrepresentable",
            ));
        }
        poll(cancelled)?;
        Ok(heat)
    }
}

#[cfg(test)]
mod tests;
