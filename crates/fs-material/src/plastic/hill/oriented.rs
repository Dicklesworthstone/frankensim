//! Fixed material axes for the same Hill return map and retained local history.
//!
//! Columns of the supplied proper rotation are the rolling, transverse and
//! thickness axes in spatial coordinates. Stress and strain are symmetric
//! tensors with tensor shear; their complete six-component transformations
//! are prepared once. Every constitutive evaluation uses the original local
//! law, including its domain, hardening, energy and analytical tangent.
//! Material history stays in these immutable local axes. This is a fixed-axis
//! small-strain model, not a finite-rotation texture evolution law.

use super::{Hill48, HillEnergyIncrement, HillError, HillResponse, HillState};
use crate::plastic::energy::{EnergyTolerances, StoredEnergy};
use crate::tensor::rotate;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Original Hill law with an explicitly chosen, immutable material orientation.
#[derive(Clone, Debug, PartialEq)]
pub struct OrientedHill48 {
    law: Hill48,
    material_to_spatial: [[f64; 3]; 3],
    to_local: Tangent6,
    to_spatial: Tangent6,
    aligned: bool,
}

impl OrientedHill48 {
    /// Bind local material axes to spatial coordinates. Reflections, scaling,
    /// shears and nonfinite rotations refuse; rotations are never normalized.
    pub fn new(law: Hill48, material_to_spatial: [[f64; 3]; 3]) -> Result<Self, HillError> {
        let q = material_to_spatial;
        if q.iter().flatten().any(|x| !x.is_finite()) {
            return Err(HillError::Parameters("Hill orientation must be finite"));
        }
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| q[k][i] * q[k][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                if !dot.is_finite() || (dot - expected).abs() > 1e-12 {
                    return Err(HillError::Parameters(
                        "Hill material axes must be orthonormal",
                    ));
                }
            }
        }
        let determinant = q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1])
            - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0])
            + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]);
        if !determinant.is_finite() || (determinant - 1.0).abs() > 1e-12 {
            return Err(HillError::Parameters(
                "Hill orientation must be a proper rotation",
            ));
        }
        let transpose = core::array::from_fn(|i| core::array::from_fn(|j| q[j][i]));
        Ok(Self {
            law,
            material_to_spatial,
            to_local: transform(&transpose),
            to_spatial: transform(&q),
            aligned: q == IDENTITY,
        })
    }

    /// The unchanged native law in material coordinates.
    pub const fn law(&self) -> &Hill48 {
        &self.law
    }

    /// Columns are the material axes in spatial coordinates.
    pub const fn material_to_spatial(&self) -> &[[f64; 3]; 3] {
        &self.material_to_spatial
    }

    /// Inspect retained plastic strain in spatial coordinates without evolving it.
    pub fn plastic_strain_spatial(&self, state: &HillState) -> Voigt {
        if self.aligned {
            state.plastic_strain
        } else {
            apply(&self.to_spatial, &state.plastic_strain)
        }
    }

    fn local_strain(&self, strain: &Voigt) -> Voigt {
        if self.aligned {
            *strain
        } else {
            apply(&self.to_local, strain)
        }
    }

    fn spatial_response(&self, mut response: HillResponse) -> Result<HillResponse, HillError> {
        if self.aligned {
            return Ok(response);
        }
        response.stress_pa = apply(&self.to_spatial, &response.stress_pa);
        // T_spatial C_local T_local preserves tensor-shear work conjugacy.
        let right: Tangent6 = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..6)
                    .map(|k| response.tangent_pa[i][k] * self.to_local[k][j])
                    .sum()
            })
        });
        response.tangent_pa = core::array::from_fn(|i| {
            core::array::from_fn(|j| (0..6).map(|k| self.to_spatial[i][k] * right[k][j]).sum())
        });
        if response
            .stress_pa
            .iter()
            .chain(response.tangent_pa.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(HillError::Numerical(
                "rotated Hill response is not representable",
            ));
        }
        Ok(response)
    }

    /// Native trial with spatial stress/tangent and unchanged local state semantics.
    pub fn response(&self, strain: &Voigt, old: &HillState) -> Result<HillResponse, HillError> {
        self.spatial_response(self.law.response(&self.local_strain(strain), old)?)
    }

    /// Inspect the original elastic/hardening storage without a second return.
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &HillState,
    ) -> Result<StoredEnergy, HillError> {
        self.law.stored_energy(&self.local_strain(strain), state)
    }

    /// Inspect the accepted local history at its actual spatial strain without
    /// advancing the return map or changing the material frame.
    pub fn validate_committed_history(
        &self,
        strain: &Voigt,
        state: &HillState,
        tolerance_pa: f64,
    ) -> Result<(), HillError> {
        self.law
            .validate_committed_history(&self.local_strain(strain), state, tolerance_pa)
    }

    /// Preserve the original work/heat account under the fixed rotation.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old: &HillState,
        strain: &Voigt,
        tolerance: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        let mut increment = self.law.advance_with_energy(
            &self.local_strain(old_strain),
            old,
            &self.local_strain(strain),
            tolerance,
        )?;
        increment.strain = *strain;
        increment.response = self.spatial_response(increment.response)?;
        Ok(increment)
    }

    pub(crate) fn advance_with_energy_from(
        &self,
        previous: &Self,
        old_strain: &Voigt,
        old: &HillState,
        strain: &Voigt,
        tolerance: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        if self.material_to_spatial != previous.material_to_spatial {
            return Err(HillError::Parameters(
                "thermal Hill history requires the original material frame",
            ));
        }
        let mut increment = self.law.advance_with_energy_from(
            &previous.law,
            &previous.local_strain(old_strain),
            old,
            &self.local_strain(strain),
            tolerance,
        )?;
        increment.strain = *strain;
        increment.response = self.spatial_response(increment.response)?;
        Ok(increment)
    }
}

fn transform(rotation: &[[f64; 3]; 3]) -> Tangent6 {
    let columns: [Voigt; 6] = core::array::from_fn(|j| {
        let mut basis = [0.0; 6];
        basis[j] = 1.0;
        rotate(&basis, rotation)
    });
    core::array::from_fn(|i| core::array::from_fn(|j| columns[j][i]))
}

fn apply(matrix: &Tangent6, value: &Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * value[j]).sum())
}

impl SmallStrainLaw for OrientedHill48 {
    type State = HillState;
    fn initial_state(&self) -> HillState {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &HillState) -> Voigt {
        // A residual evaluation needs only the returned stress. Avoid forming
        // and rotating the full constitutive tangent during each line search.
        let local = self.law.stress(&self.local_strain(strain), state);
        let stress = if self.aligned {
            local
        } else {
            apply(&self.to_spatial, &local)
        };
        if stress.iter().all(|v| v.is_finite()) {
            stress
        } else {
            [f64::NAN; 6]
        }
    }
    fn tangent(&self, strain: &Voigt, state: &HillState) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |r| r.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &HillState) -> HillState {
        // Accepted history remains local; neither stress nor tangent needs a
        // spatial rotation for this operation.
        self.law.update_state(&self.local_strain(strain), state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}
