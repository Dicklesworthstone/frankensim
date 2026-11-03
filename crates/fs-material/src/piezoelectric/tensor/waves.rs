//! Local piezoacoustic waves from the same material used by static transducers.
//!
//! In a charge-free insulating bulk, electroquasistatics requires E parallel to
//! the wave normal n and n.D=0. Eliminating that ONE electric coordinate gives
//! `C_eff(n)=C^E+(e^T n)(e^T n)^T/(n.epsilon^S.n)`. In particular this is NOT
//! the full fixed-D tensor C^D: transverse D need not vanish. The existing
//! elastic Christoffel/eigenmode owner evaluates this direction-dependent law.
//! See COMSOL's *Piezoelectricity*, Wave Propagation in Piezoelectric Media:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.043.html>.
//!
//! Phase speeds and energy velocities belong to the frozen lossless model.
//! For charge-free waves the energy velocity differentiates the electrical
//! constraint as the wave normal changes; freezing C_eff(n) misses electric
//! energy transport transverse to the wave normal.
//! These are infinitesimal local wave properties,
//! not a finite-device resonance solve, a frequency-dependent material fit,
//! acoustic radiation or a measured ultrasonic calibration.

use super::{Result, TensorPiezoelectric, invalid};
use crate::elastic::acoustic::{DirectionalElasticWaves, ElasticAcousticState};
use crate::tensor::mandel_basis;

/// Electrical constraint of the local plane-wave calculation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PiezoelectricWaveCondition {
    /// Electric-field perturbation is clamped to zero; uses C^E.
    FixedField,
    /// No bulk free charge, curl(E)=0: only n.D vanishes.
    ChargeFree,
}

/// The native mechanical eigenmodes and their reciprocal electric response.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectionalPiezoelectricWaves {
    /// Explicit electrical condition; device electrode constraints are separate.
    pub condition: PiezoelectricWaveCondition,
    /// Three native Christoffel modes and source identity, in ascending speed
    /// order. Energy velocity includes the reciprocal electric contribution
    /// for charge-free waves, not only mechanical stress power.
    pub mechanical: DirectionalElasticWaves,
    /// Electric field per unit displacement-gradient amplitude [V/m], one
    /// vector per native mode. For a plane wave of displacement amplitude A
    /// and wavenumber k, multiply by k*A with the strain's phase convention.
    pub field_per_displacement_gradient_v_m: [[f64; 3]; 3],
}

impl TensorPiezoelectric {
    /// Compute local acoustic modes with the SAME density, stiffness, coupling
    /// and dielectric tensor as the mechanical/electric energy. The explicit
    /// strain ceiling is passed to the existing native acoustic material owner;
    /// no finite-amplitude or frequency-domain validity is added by this query.
    pub fn directional_waves(
        &self,
        direction: [f64; 3],
        condition: PiezoelectricWaveCondition,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<DirectionalPiezoelectricWaves> {
        let scale = direction.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        if direction.iter().any(|x| !x.is_finite())
            || scale == 0.0
            || !maximum_absolute_mandel_strain.is_finite()
            || maximum_absolute_mandel_strain <= 0.0
            || maximum_absolute_mandel_strain >= 1.0
        {
            return Err(invalid(
                "finite nonzero wave direction and small-signal strain ceiling in (0,1) required",
            ));
        }
        let n = direction.map(|x| x / scale);
        let length = fs_math::det::sqrt(n.iter().map(|x| x * x).sum());
        let n = n.map(|x| x / length);
        let p: [f64; 6] =
            core::array::from_fn(|a| (0..3).map(|i| n[i] * self.coupling_c_m2()[i][a]).sum());
        let epsilon = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| n[i] * self.permittivity_f_m()[i][j] * n[j])
                    .sum::<f64>()
            })
            .sum::<f64>();
        if !epsilon.is_finite() || epsilon <= 0.0 {
            return Err(invalid(
                "directional dielectric response is not positive and representable",
            ));
        }
        let mut effective = *self.stiffness_mandel_pa();
        if condition == PiezoelectricWaveCondition::ChargeFree {
            for i in 0..6 {
                for j in i..6 {
                    let value = effective[i][j] + p[i] * (p[j] / epsilon);
                    if !value.is_finite() {
                        return Err(invalid("piezoacoustic stiffening is not representable"));
                    }
                    effective[i][j] = value;
                    effective[j][i] = value;
                }
            }
        }
        let acoustic = ElasticAcousticState::try_new_mandel(
            effective,
            self.density_kg_m3(),
            self.identity(),
            maximum_absolute_mandel_strain,
        )
        .map_err(|_| invalid("native acoustic material admission failed"))?;
        let mut mechanical = acoustic
            .directional_waves(n)
            .map_err(|_| invalid("native piezoacoustic eigenmode solve failed"))?;
        let mut fields = [[0.0; 3]; 3];
        if condition == PiezoelectricWaveCondition::ChargeFree {
            let basis = mandel_basis();
            for (mode, field) in fields.iter_mut().enumerate() {
                let displacement = mechanical.modes[mode].polarization;
                let strain: [f64; 6] = core::array::from_fn(|a| {
                    (0..3)
                        .map(|i| {
                            (0..3)
                                .map(|j| basis[a][i][j] * n[j] * displacement[i])
                                .sum::<f64>()
                        })
                        .sum()
                });
                let scalar = -(0..6).map(|a| p[a] * strain[a]).sum::<f64>() / epsilon;
                *field = n.map(|x| scalar * x);
                if field.iter().any(|x| !x.is_finite()) {
                    return Err(invalid("acoustic electric field is not representable"));
                }
                // Write t=(p.strain)/epsilon=-E_n. Differentiating the
                // homogeneous wave-vector energy adds t*D/(rho*c) to the
                // frozen-C_eff mechanical velocity, where
                // D=e:strain-t*epsilon*n and n.D=0. Both the coupling and
                // dielectric normal must vary; a fixed-D tensor is different.
                let displacement: [f64; 3] = core::array::from_fn(|i| {
                    (0..6)
                        .map(|a| self.coupling_c_m2()[i][a] * strain[a])
                        .sum::<f64>()
                        + (0..3)
                            .map(|j| self.permittivity_f_m()[i][j] * field[j])
                            .sum::<f64>()
                });
                let mode = &mut mechanical.modes[mode];
                for (velocity, electric_displacement) in
                    mode.energy_velocity_m_s.iter_mut().zip(displacement)
                {
                    *velocity -= (scalar / mode.impedance_pa_s_m) * electric_displacement;
                }
                if mode.energy_velocity_m_s.iter().any(|v| !v.is_finite()) {
                    return Err(invalid(
                        "piezoacoustic energy velocity is not representable",
                    ));
                }
            }
        }
        Ok(DirectionalPiezoelectricWaves {
            condition,
            mechanical,
            field_per_displacement_gradient_v_m: fields,
        })
    }
}
