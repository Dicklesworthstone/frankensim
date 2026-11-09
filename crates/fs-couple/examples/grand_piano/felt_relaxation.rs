//! Thin mechanical adapter for fs-material's original passive creep owner.
//!
//! Prony realization, branch propagation and independent heat belong to
//! `fiber::felt_rate::CreepSpectrum`. This adapter supplies only the contact
//! area's and thickness's work-conjugate SI map: stress=force/area,
//! strain=deformation/thickness, and energy=energy_density*area*thickness.
//! The piano's original WoolFelt discrete-gradient force solve is unchanged.
//! Its held force drives these SAME histories; memory recovers while separated.
//! No allocation or branch exponential occurs in the warm adapter methods.
use fs_material::fiber::felt_rate::{CreepSpectrum, MAX_FELT_RATE_BRANCHES, PreparedCreep};
use fs_material::visco::GeneralizedMaxwell;

pub const MAX_TERMS: usize = MAX_FELT_RATE_BRANCHES;

/// Mechanical creep coordinates [m], retained for existing checkpoint consumers.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Memory(pub [f64; MAX_TERMS]);

#[derive(Clone, Debug)]
pub struct Prepared {
    material: PreparedCreep,
    area_m2: f64,
    thickness_m: f64,
    volume_m3: f64,
    compliance_m_n: f64,
}

#[derive(Clone, Debug)]
pub struct Spectrum {
    material: CreepSpectrum,
}

/// Explicitly estimated felt time scales. Instantaneous reference modulus is
/// 5 MPa, matching the virgin tangent of demonstration_law at strain 0.2.
/// These are NOT a measured Steinway hammer coupon or a claimed Stulov fit.
pub fn demonstration_prony() -> GeneralizedMaxwell {
    GeneralizedMaxwell::new(2.5e6, vec![(2.0e6, 0.0002), (0.5e6, 0.004)])
        .expect("admissible authored Prony constants")
}

impl Spectrum {
    /// Use the shared material source conversion, with no fitted envelope or
    /// dropped weak weights. Unresolved original spectra refuse at this owner.
    pub fn from_prony(card: &GeneralizedMaxwell) -> Result<Self, String> {
        CreepSpectrum::from_prony(card)
            .map(|material| Self { material })
            .map_err(|error| error.to_string())
    }

    pub fn prepare(&self, area_m2: f64, thickness_m: f64, dt: f64) -> Result<Prepared, String> {
        if [area_m2, thickness_m, dt]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err("creep preparation needs positive SI area/thickness/dt".into());
        }
        let material = self
            .material
            .prepare(dt)
            .map_err(|error| error.to_string())?;
        let volume_m3 = area_m2 * thickness_m;
        let deformation_per_stress = thickness_m / area_m2;
        let compliance_m_n = material.compliance_pa_inverse() * deformation_per_stress;
        if !volume_m3.is_finite()
            || volume_m3 <= 0.0
            || !deformation_per_stress.is_finite()
            || deformation_per_stress <= 0.0
            || !compliance_m_n.is_finite()
            || (material.compliance_pa_inverse() > 0.0 && compliance_m_n <= 0.0)
        {
            return Err("creep contact geometry or compliance is unrepresentable".into());
        }
        Ok(Prepared {
            material,
            area_m2,
            thickness_m,
            volume_m3,
            compliance_m_n,
        })
    }
}

impl Prepared {
    pub fn compliance(&self) -> f64 {
        self.compliance_m_n
    }

    fn strains(&self, memory: &Memory) -> [f64; MAX_TERMS] {
        memory.0.map(|z| z / self.thickness_m)
    }

    /// Deformation of an admitted prepared history [m].
    ///
    /// Existing live consumers own and retain admitted histories. A malformed
    /// manually constructed checkpoint is a programming error for this legacy
    /// query; the fallible update below rejects such input before publication.
    pub fn deformation(&self, memory: &Memory) -> f64 {
        self.material
            .strain(&self.strains(memory))
            .expect("admitted creep history")
            * self.thickness_m
    }

    /// Zero-force predicted deformation used by the SAME implicit contact solve.
    pub fn free_deformation(&self, memory: &Memory) -> f64 {
        self.material
            .free_strain(&self.strains(memory))
            .expect("admitted creep history")
            * self.thickness_m
    }

    pub fn stored(&self, memory: &Memory) -> f64 {
        self.material
            .stored_energy_density(&self.strains(memory))
            .expect("admitted creep history")
            * self.volume_m3
    }

    /// Trial only: the original material owner advances the common held force.
    /// The input is never modified; commit only after the coupled step passes.
    pub fn try_advance(&self, memory: &Memory, force: f64) -> Result<(Memory, f64), &'static str> {
        if !force.is_finite() || force < 0.0 {
            return Err("felt creep needs finite nonnegative contact force");
        }
        let result = self
            .material
            .advance(&self.strains(memory), force / self.area_m2)
            .map_err(|_| "shared felt creep refused the original material history or work")?;
        let next = Memory(result.strain.map(|strain| strain * self.thickness_m));
        let heat = result.dissipated_energy_j_m3 * self.volume_m3;
        if next.0.iter().any(|x| !x.is_finite()) || !heat.is_finite() || heat < 0.0 {
            return Err("felt creep contact state or heat is unrepresentable");
        }
        Ok((next, heat))
    }

    /// Compatibility entry point for existing live consumers with admitted
    /// histories and finite nonnegative solved forces. Reference owners use
    /// `try_advance` to propagate numerical refusal through their transaction.
    ///
    /// # Panics
    /// If a manually supplied history, force or numerical response is refused.
    pub fn advance(&self, memory: &Memory, force: f64) -> (Memory, f64) {
        self.try_advance(memory, force)
            .expect("admitted prepared creep trial")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-200),
            "{a:e} != {b:e} (relative tolerance {relative:e})"
        );
    }

    #[test]
    fn creep_realization_matches_its_existing_prony_owner() {
        let gm = demonstration_prony();
        let spectrum = Spectrum::from_prony(&gm).unwrap();
        let instant = gm.e_inf + gm.terms.iter().map(|t| t.0).sum::<f64>();
        for i in 0..=80 {
            let w = 10.0f64.powf(-1.0 + 8.0 * i as f64 / 80.0);
            let (mut jr, mut ji) = (1.0 / instant, 0.0);
            for branch in spectrum.material.branches() {
                let x = w * (branch.viscosity_pa_s / branch.modulus_pa);
                let response = 1.0 / (branch.modulus_pa * (1.0 + x * x));
                jr += response;
                ji -= x * response;
            }
            let (er, ei) = gm.modulus(w);
            let den = er * er + ei * ei;
            close(jr, er / den, 1e-12);
            close(ji, -ei / den, 1e-12);
        }
    }

    #[test]
    fn held_force_and_release_close_work_without_negative_loss() {
        let spectrum = Spectrum::from_prony(&demonstration_prony()).unwrap();
        for dt in [1e-8, 1.0 / 192_000.0, 1.0 / 8_000.0, 0.1] {
            let prepared = spectrum.prepare(1e-4, 0.008, dt).unwrap();
            let mut memory = Memory::default();
            let (mut work, mut dissipated) = (0.0, 0.0);
            for i in 0..2000 {
                let force = if i < 1000 {
                    30.0 * (i as f64 * 0.03).sin().abs()
                } else {
                    0.0
                };
                let (next, loss) = prepared.try_advance(&memory, force).unwrap();
                work += force * (prepared.deformation(&next) - prepared.deformation(&memory));
                dissipated += loss;
                assert!(loss >= 0.0);
                assert!(next.0.iter().all(|z| *z >= 0.0));
                assert!((work - prepared.stored(&next) - dissipated).abs() < 1e-9);
                memory = next;
            }
        }
    }

    #[test]
    fn analytic_material_creep_and_recovery_have_work_conjugate_geometry() {
        let card = GeneralizedMaxwell::new(1e6, vec![(2e6, 0.002)]).unwrap();
        let spectrum = Spectrum::from_prony(&card).unwrap();
        let (modulus, tau) = (1.5e6, 0.006);
        let (dt, stress) = (0.001, 2000.0);
        for (area, thickness) in [(1e-4, 0.008), (2e-4, 0.008), (1e-4, 0.016)] {
            let prepared = spectrum.prepare(area, thickness, dt).unwrap();
            let force = stress * area;
            let k = modulus * area / thickness;
            let mut memory = Memory::default();
            let mut heat = 0.0;
            for step in 1..=20 {
                let (next, loss) = prepared.try_advance(&memory, force).unwrap();
                heat += loss;
                memory = next;
                let decay = (-dt * step as f64 / tau).exp();
                close(
                    prepared.deformation(&memory),
                    force / k * (1.0 - decay),
                    1e-12,
                );
                close(heat, 0.5 * force * force / k * (1.0 - decay * decay), 1e-12);
            }
            let before = prepared.deformation(&memory);
            let (next, loss) = prepared.try_advance(&memory, 0.0).unwrap();
            let decay = (-dt / tau).exp();
            close(prepared.deformation(&next), before * decay, 1e-12);
            close(
                loss,
                0.5 * k * before * before * (1.0 - decay * decay),
                1e-12,
            );
            close(
                prepared.stored(&next),
                0.5 * k * prepared.deformation(&next).powi(2),
                1e-12,
            );
        }
    }

    #[test]
    fn equal_stress_area_partition_preserves_displacement_and_total_heat() {
        let spectrum = Spectrum::from_prony(&demonstration_prony()).unwrap();
        let whole = spectrum.prepare(1e-4, 0.008, 1.0 / 192_000.0).unwrap();
        let fractions = [0.1, 0.2, 0.3, 0.4];
        let sites = fractions.map(|fraction| {
            spectrum
                .prepare(fraction * 1e-4, 0.008, 1.0 / 192_000.0)
                .unwrap()
        });
        let mut original = Memory::default();
        let mut memories = [Memory::default(); 4];
        for i in 0..500 {
            let force = if i < 300 {
                25.0 * (i as f64 * 0.03).sin().abs()
            } else {
                0.0
            };
            let (next, heat) = whole.try_advance(&original, force).unwrap();
            let mut sum_heat = 0.0;
            let mut sum_storage = 0.0;
            for j in 0..4 {
                let (local, lost) = sites[j]
                    .try_advance(&memories[j], fractions[j] * force)
                    .unwrap();
                close(
                    sites[j].deformation(&local),
                    whole.deformation(&next),
                    2e-12,
                );
                sum_heat += lost;
                sum_storage += sites[j].stored(&local);
                memories[j] = local;
            }
            close(sum_heat, heat, 2e-12);
            close(sum_storage, whole.stored(&next), 2e-12);
            original = next;
        }
    }

    #[test]
    fn positive_weak_source_loss_and_long_recovery_tail_are_retained() {
        // This source strength was below the old adapter's weight*E0 cutoff.
        let card = GeneralizedMaxwell::new(1e6, vec![(1e-12, 0.2)]).unwrap();
        let spectrum = Spectrum::from_prony(&card).unwrap();
        let prepared = spectrum.prepare(1e-4, 0.008, 0.1).unwrap();
        let (next, heat) = prepared.try_advance(&Memory::default(), 3.0).unwrap();
        assert!(prepared.compliance() > 0.0 && prepared.deformation(&next) > 0.0 && heat > 0.0);
        let k = 1e24 * 1e-4 / 0.008;
        let decay = (-0.5_f64).exp();
        close(prepared.deformation(&next), 3.0 / k * (1.0 - decay), 1e-12);
        close(heat, 0.5 * 9.0 / k * (1.0 - decay * decay), 1e-12);

        let card = GeneralizedMaxwell::new(1e6, vec![(2e6, 0.002)]).unwrap();
        let prepared = Spectrum::from_prony(&card)
            .unwrap()
            .prepare(1e-4, 0.008, 0.6)
            .unwrap();
        let (loaded, _) = prepared.try_advance(&Memory::default(), 3.0).unwrap();
        let (recovered, _) = prepared.try_advance(&loaded, 0.0).unwrap();
        assert!(prepared.deformation(&recovered) > 0.0);
        close(
            prepared.deformation(&recovered),
            prepared.deformation(&loaded) * (-100.0_f64).exp(),
            1e-12,
        );
    }

    #[test]
    fn elastic_card_has_no_relaxation_memory_and_zero_terms_are_legal() {
        let card = GeneralizedMaxwell::new(1e6, vec![(0.0, 0.1)]).unwrap();
        let prepared = Spectrum::from_prony(&card)
            .unwrap()
            .prepare(1.0, 1.0, 0.01)
            .unwrap();
        assert_eq!(prepared.compliance(), 0.0);
        assert_eq!(
            prepared.deformation(&Memory::default()).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            prepared.free_deformation(&Memory::default()).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            prepared.stored(&Memory::default()).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            prepared.try_advance(&Memory::default(), 3.0).unwrap(),
            (Memory::default(), 0.0)
        );
        assert_eq!(
            prepared.advance(&Memory::default(), 3.0),
            (Memory::default(), 0.0)
        );
        assert!(
            Spectrum::from_prony(&GeneralizedMaxwell {
                e_inf: f64::NAN,
                terms: vec![]
            })
            .is_err()
        );
    }

    #[test]
    fn refused_material_trials_preserve_mechanical_memory_and_retry() {
        let prepared = Spectrum::from_prony(&demonstration_prony())
            .unwrap()
            .prepare(1e-4, 0.008, 1.0 / 192_000.0)
            .unwrap();
        let initial = Memory::default();
        let expected = prepared.try_advance(&initial, 3.0).unwrap();
        assert!(prepared.try_advance(&initial, -1.0).is_err());
        assert!(prepared.try_advance(&initial, f64::NAN).is_err());
        let mut invalid = initial;
        invalid.0[7] = 0.1;
        assert!(prepared.try_advance(&invalid, 3.0).is_err());
        assert_eq!(initial, Memory::default());
        assert_eq!(prepared.try_advance(&initial, 3.0).unwrap(), expected);
    }
}
