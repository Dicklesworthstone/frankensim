//! Source-bound, fully relaxed initial conditions for original material histories.
//!
//! Initial storage is supplied as an initial condition. No loading trajectory,
//! prestress solve, elapsed relaxation time, work input or heat is invented.

use super::{Result, TensorPronyHistory, Voigt, invalid, mandel, poll};
use crate::plane_stress_viscosity::MaterialPlane;
use crate::visco::tensor::HarmonicWorkSum;
use fs_math::{c64::C64, dd::Dd};

fn finite(value: Dd) -> Result<Dd> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(invalid("relaxed initial material state is unrepresentable"))
    }
}

impl TensorPronyHistory {
    fn require_virgin_template(&self, cancelled: &mut impl FnMut() -> bool) -> Result<()> {
        poll(cancelled)?;
        if self.strain.iter().chain(&self.stress).any(|v| *v != 0.0)
            || self.stored != 0.0
            || self.work != 0.0
            || self.heat != 0.0
        {
            return Err(invalid("relaxed initialization requires an original virgin history template"));
        }
        for memory in &self.q {
            poll(cancelled)?;
            if memory.iter().any(|v| *v != 0.0) {
                return Err(invalid("relaxed initialization cannot discard original Maxwell memory"));
            }
        }
        Ok(())
    }

    fn relaxed_mandel_initial(
        &self,
        strain: Voigt,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        self.model.validate_mandel_strain(&strain)
            .map_err(|_| invalid("relaxed initial strain leaves the original material domain"))?;
        let mut elastic = [Dd::ZERO; 6];
        for (i, row) in self.model.equilibrium_gram.iter().enumerate() {
            let mut sum = HarmonicWorkSum::default();
            for j in 0..6 {
                sum.add([row[j], strain[j], 1.0, 1.0]);
            }
            elastic[i] = finite(sum.retained_value())?;
        }
        let mut stress = [0.0; 6];
        let mut storage = HarmonicWorkSum::default();
        for i in 0..6 {
            let mut sum = HarmonicWorkSum::default();
            for j in 0..6 {
                for word in [elastic[j].hi, elastic[j].lo] {
                    sum.add([self.model.equilibrium_gram[j][i], word, 1.0, 1.0]);
                }
            }
            stress[i] = finite(sum.retained_value())?.to_f64();
            for a in [elastic[i].hi, elastic[i].lo] {
                for b in [elastic[i].hi, elastic[i].lo] {
                    storage.add([0.5, a, b, 1.0]);
                }
            }
        }
        let stored = finite(storage.retained_value())?.to_f64();
        if !stored.is_finite() || stored < 0.0 || stress.iter().any(|v| !v.is_finite()) {
            return Err(invalid("relaxed initial stress or storage is unrepresentable"));
        }
        poll(cancelled)?;
        // The original Arc and all branch rows survive. The checked virgin
        // template already has q=0 and zero cumulative work/heat.
        let mut next = self.clone();
        next.strain = strain;
        next.stress = stress;
        next.stored = stored;
        poll(cancelled)?;
        Ok(next)
    }

    /// Create an explicitly fully relaxed initial condition at the supplied
    /// tensor-shear Voigt strain, sharing this virgin template's exact source.
    ///
    /// All original Maxwell branch coordinates are zero. Initial strain rate is
    /// zero, so Kelvin viscosity adds no initial stress or heat. Initial storage
    /// is the original equilibrium spring energy; cumulative work and heat start
    /// at zero. The caller supplies this initial energy and any balancing load.
    /// No preceding loading path, time to relax, or mechanical equilibrium of an
    /// assembled body is inferred. This is initialization, not a history update:
    /// nonvirgin templates refuse, and this borrowed template is never changed.
    ///
    /// Nonzero static strain requires the original source to admit DC. A zero
    /// initial strain remains valid for an AC-only law. Original full-strain
    /// admission, representability and cancellation apply before publication.
    pub fn initial_relaxed_strain(
        &self,
        strain: Voigt,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        self.require_virgin_template(&mut cancelled)?;
        let strain = mandel(&strain);
        if strain.iter().any(|v| !v.is_finite()) {
            return Err(invalid("relaxed initial strain must be finite"));
        }
        if strain.iter().any(|v| *v != 0.0) {
            self.model.modulus_mandel(0.0)?;
        }
        self.relaxed_mandel_initial(strain, &mut cancelled)
    }

    /// Create a fully relaxed plane-stress initial condition on the SAME source.
    ///
    /// Input is ordered engineering [epsilon_aa, epsilon_bb, gamma_ab] in
    /// `plane`. The original DC plane-stress owner recovers every transverse
    /// strain from the complete equilibrium tensor. No branch is condensed,
    /// populated or discarded, and no finite-duration evolution is performed.
    /// The returned full strain/stress use tensor shears through their getters.
    ///
    /// The source, initial-energy and virgin-template rules of
    /// [`Self::initial_relaxed_strain`] apply. This permits a restrained thermal
    /// preload when the caller subtracts its independently admitted held free
    /// strain from the compatible total in-plane strain.
    pub fn initial_relaxed_plane_stress(
        &self,
        plane: MaterialPlane,
        in_plane_engineering: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        self.require_virgin_template(&mut cancelled)?;
        if in_plane_engineering.iter().any(|v| !v.is_finite()) {
            return Err(invalid("relaxed initial plane strain must be finite"));
        }
        if in_plane_engineering.iter().all(|v| *v == 0.0) {
            poll(&mut cancelled)?;
            return Ok(self.clone());
        }
        let recovered = self.model.prepare_plane_stress_harmonic(plane, 0.0)
            .and_then(|prepared| prepared.recover_strain(&in_plane_engineering.map(C64::from_re)))
            .map_err(|_| invalid("relaxed plane-stress initialization leaves the original DC or strain domain"))?;
        if recovered.iter().any(|v| !v.re.is_finite() || v.im != 0.0) {
            return Err(invalid("relaxed plane-stress initial strain is not real and finite"));
        }
        let mut next = self.relaxed_mandel_initial(recovered.map(|v| v.re), &mut cancelled)?;
        let outside = match plane {
            MaterialPlane::Axes12 | MaterialPlane::Axes21 => [2, 4, 5],
            MaterialPlane::Axes13 => [1, 3, 4],
            MaterialPlane::Axes23 => [0, 3, 5],
        };
        for i in outside {
            let mut scale = HarmonicWorkSum::default();
            for j in 0..6 {
                scale.add([self.model.c_eq[i][j].abs(), next.strain[j].abs(), 1.0, 1.0]);
            }
            let scale = finite(scale.retained_value())?.to_f64();
            if !scale.is_finite() || next.stress[i].abs() > 4096.0 * f64::EPSILON * scale {
                return Err(invalid("original relaxed transverse traction is unresolved"));
            }
            next.stress[i] = 0.0;
        }
        poll(&mut cancelled)?;
        Ok(next)
    }
}
