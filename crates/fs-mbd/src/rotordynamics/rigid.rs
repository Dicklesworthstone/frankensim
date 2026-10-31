//! A nominal axisymmetric rotor on finite, possibly cross-coupled radial bearings.
//!
//! Coordinates `[x, y, theta_x, theta_y]` describe small lateral translations
//! and right-handed tilts in STATIONARY reference axes, with positive spin
//! about +z. The reference plane passes through the nominal mass centroid.
//! The original full body inertia is lowered into these axes once. A bearing
//! at z sees `[x + z theta_y, y - z theta_x]`; its generalized force is the
//! exact transpose of that displacement map, including bending moments.
//!
//! For `n = [theta_y, -theta_x, 1]`, transverse angular momentum is
//! `H = Jd theta_dot + Jp Omega n`. Its derivative gives
//! `G[theta_x,theta_y] = +Jp` and `G[theta_y,theta_x] = -Jp`.
//! No gyroscopic term is assigned heat or symmetrized away. Constant-speed
//! unbalance is a declared small perturbation of the nominal rotor, not an
//! additional mass/inertia insertion. The constant-speed Jeffcott convention
//! is described in NPTEL, Simple Rotor Systems, lecture 2:
//! <https://archive.nptel.ac.in/content/storage2/courses/112103024/module2/lec2/2.html>.
//!
//! Bearings are actual coefficients at one operating point. Their stiffness
//! and damping may be anisotropic, nonsymmetric or active. Changing speed,
//! temperature, preload or a fluid state requires supplying the corresponding
//! coefficients. This is a rigid-rotor linearization; shaft flexibility,
//! nonlinear bearing clearance, spin acceleration and asymmetric-rotor Floquet
//! physics require their original dynamic operators, not this reduction.

use super::{RotorBearingModel, RotorError};
use crate::articulated::SpatialInertia;
use fs_ga::{So3, Vec3};
use fs_math::c64::C64;

/// Bound on supplied bearing and unbalance stations per rigid rotor query.
pub const MAX_ROTOR_STATIONS: usize = 256;

/// One linear radial bearing about its actual operating equilibrium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialBearing {
    /// Signed axial distance from the nominal centroid plane [m].
    pub station_m: f64,
    /// Full `[x,y]` stiffness in the stationary reference axes [N/m].
    pub stiffness_n_per_m: [[f64; 2]; 2],
    /// Full `[x,y]` velocity-force coefficients in the same axes [N s/m].
    pub damping_ns_per_m: [[f64; 2]; 2],
}

/// A prescribed residual mass-eccentricity vector rotating at the shaft speed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorUnbalance {
    /// Signed axial distance from the nominal centroid plane [m].
    pub station_m: f64,
    /// Nonnegative mass times eccentric radius [kg m].
    pub mass_eccentricity_kg_m: f64,
    /// Direction from +x towards +y at physical time zero [rad].
    pub phase_rad: f64,
}

/// Actual complex bearing displacement and force on the rotor.
#[derive(Clone, Debug, PartialEq)]
pub struct BearingReaction {
    /// Original supplied bearing index.
    pub bearing: usize,
    /// Radial displacement at the bearing [m].
    pub displacement_m: [C64; 2],
    /// Radial force exerted BY the bearing ON the rotor [N].
    pub force_n: [C64; 2],
    /// Bending moment about the centroid plane, `[tau_x,tau_y]` [N m].
    pub moment_nm: [C64; 2],
    /// Signed cycle-average bearing work on the rotor [W]. This is not an
    /// assigned heat rate for arbitrary nonsymmetric or active coefficients.
    pub mean_power_on_rotor_w: f64,
}

/// Source inertia and physical stations retained with their assembled operators.
#[derive(Clone, Debug)]
pub struct RigidRotor {
    inertia: SpatialInertia,
    body_from_rotor: So3,
    axisymmetry_relative_residual: f64,
    bearings: Vec<RadialBearing>,
    model: RotorBearingModel,
}

impl RigidRotor {
    /// Lower the original centroidal tensor into a supplied proper reference
    /// frame and assemble all bearing force/moment maps by virtual work.
    ///
    /// `body_from_rotor` maps reference rotor axes into original material axes
    /// at the nominal configuration; the perturbation axes do not co-rotate.
    /// Axisymmetry must hold to the explicit relative tolerance in `(0,1e-8]`.
    /// The transverse moment uses the mean of the two admitted equal moments;
    /// the measured tensor discrepancy remains queryable. Asymmetric inertia
    /// refuses instead of silently replacing its time-periodic dynamics.
    ///
    /// Zero bearings are allowed: they retain the physical free modes. The
    /// native harmonic solver refuses singular response; no support floor is
    /// inserted. Inputs remain unchanged on cancellation or invalid geometry.
    pub fn new(
        inertia: SpatialInertia,
        body_from_rotor: So3,
        bearings: &[RadialBearing],
        axisymmetry_relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RotorError> {
        poll(&mut cancelled)?;
        if bearings.len() > MAX_ROTOR_STATIONS
            || !axisymmetry_relative_tolerance.is_finite()
            || axisymmetry_relative_tolerance <= 0.0
            || axisymmetry_relative_tolerance > 1e-8
        {
            return Err(RotorError::Invalid(
                "rigid rotor station count or axisymmetry tolerance",
            ));
        }
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
        .map(|axis| body_from_rotor.rotate(axis))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| RotorError::Invalid("rigid rotor reference frame"))?;
        let tensor = inertia.inertia_com();
        let mut lowered = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                lowered[i][j] = finite(axes[i].dot(tensor.apply(axes[j])))?;
            }
        }
        let transverse = finite(0.5 * lowered[0][0] + 0.5 * lowered[1][1])?;
        let polar = lowered[2][2];
        let scale = lowered
            .iter()
            .flatten()
            .fold(0.0_f64, |s, x| s.max(x.abs()));
        let mut defect = 0.0_f64;
        for (i, row) in lowered.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let expected = if i != j {
                    0.0
                } else if i == 2 {
                    polar
                } else {
                    transverse
                };
                defect = defect.max((value - expected).abs());
            }
        }
        let residual = finite(defect / scale)?;
        if transverse <= 0.0 || polar <= 0.0 || residual > axisymmetry_relative_tolerance {
            return Err(RotorError::Invalid(
                "source inertia is not axisymmetric about the declared spin axis",
            ));
        }
        let mut mass = [0.0; 16];
        for (i, value) in [inertia.mass(), inertia.mass(), transverse, transverse]
            .into_iter()
            .enumerate()
        {
            mass[5 * i] = value;
        }
        let mut gyro = [0.0; 16];
        gyro[11] = polar;
        gyro[14] = -polar;
        let mut stiffness = [0.0; 16];
        let mut damping = [0.0; 16];
        let mut k_correction = [0.0; 16];
        let mut c_correction = [0.0; 16];
        for bearing in bearings {
            poll(&mut cancelled)?;
            finite(bearing.station_m)?;
            let map = point_map(bearing.station_m);
            for i in 0..4 {
                for j in 0..4 {
                    for a in 0..2 {
                        for b in 0..2 {
                            let k = finite(
                                map[a][i] * finite(bearing.stiffness_n_per_m[a][b])? * map[b][j],
                            )?;
                            let c = finite(
                                map[a][i] * finite(bearing.damping_ns_per_m[a][b])? * map[b][j],
                            )?;
                            add_compensated(
                                &mut stiffness[4 * i + j],
                                &mut k_correction[4 * i + j],
                                k,
                            )?;
                            add_compensated(
                                &mut damping[4 * i + j],
                                &mut c_correction[4 * i + j],
                                c,
                            )?;
                        }
                    }
                }
            }
        }
        for i in 0..16 {
            stiffness[i] = finite(stiffness[i] + k_correction[i])?;
            damping[i] = finite(damping[i] + c_correction[i])?;
        }
        let model = RotorBearingModel::new(4, &mass, &damping, &stiffness, &gyro, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            inertia,
            body_from_rotor,
            axisymmetry_relative_residual: residual,
            bearings: bearings.to_vec(),
            model,
        })
    }

    /// Original complete mass/centroid/tensor record, without replacing its axes.
    pub const fn inertia(&self) -> SpatialInertia {
        self.inertia
    }
    /// Proper reference-to-original-material frame used for tensor lowering.
    pub const fn body_from_rotor(&self) -> So3 {
        self.body_from_rotor
    }
    /// Relative maximum-entry discrepancy in the admitted axisymmetric reduction.
    pub const fn axisymmetry_relative_residual(&self) -> f64 {
        self.axisymmetry_relative_residual
    }
    /// Full assembled model for poles, Campbell samples and forced response.
    pub fn model(&self) -> &RotorBearingModel {
        &self.model
    }

    /// Synchronous generalized force for `exp(+i Omega t)`, in `[N,N,Nm,Nm]`.
    /// Use the SAME signed speed as the harmonic forcing frequency. At phase
    /// zero the force points +x and rotates towards +y when speed is positive.
    /// Each station contributes its actual transverse force and bending moment.
    /// The supplied unbalance is a forcing perturbation; this function never
    /// inserts another mass or alters the retained nominal inertia.
    pub fn unbalance_force(
        &self,
        spin_rad_s: f64,
        unbalances: &[RotorUnbalance],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[C64; 4], RotorError> {
        poll(&mut cancelled)?;
        finite(spin_rad_s)?;
        if unbalances.len() > MAX_ROTOR_STATIONS {
            return Err(RotorError::Invalid("rotor unbalance station count"));
        }
        let mut sums = [0.0; 8];
        let mut correction = [0.0; 8];
        for point in unbalances {
            poll(&mut cancelled)?;
            finite(point.station_m)?;
            finite(point.phase_rad)?;
            if !point.mass_eccentricity_kg_m.is_finite() || point.mass_eccentricity_kg_m < 0.0 {
                return Err(RotorError::Invalid(
                    "nonnegative finite rotor mass-eccentricity required",
                ));
            }
            let amplitude = finite((point.mass_eccentricity_kg_m * spin_rad_s) * spin_rad_s)?;
            if amplitude == 0.0 && point.mass_eccentricity_kg_m != 0.0 && spin_rad_s != 0.0 {
                return Err(RotorError::Invalid(
                    "unrepresentable rotor unbalance amplitude",
                ));
            }
            let sin = fs_math::det::sin(point.phase_rad);
            let cos = fs_math::det::cos(point.phase_rad);
            let x = C64::new(amplitude * cos, amplitude * sin);
            let y = C64::new(x.im, -x.re);
            let values = [x, y, y.scale(-point.station_m), x.scale(point.station_m)];
            for (i, value) in values.into_iter().enumerate() {
                add_compensated(&mut sums[2 * i], &mut correction[2 * i], finite(value.re)?)?;
                add_compensated(
                    &mut sums[2 * i + 1],
                    &mut correction[2 * i + 1],
                    finite(value.im)?,
                )?;
            }
        }
        for i in 0..8 {
            sums[i] = finite(sums[i] + correction[i])?;
        }
        poll(&mut cancelled)?;
        Ok(std::array::from_fn(|i| {
            C64::new(sums[2 * i], sums[2 * i + 1])
        }))
    }

    /// Recover bearing forces and their physical bending moments from the
    /// accepted harmonic state. Velocity is the actual phasor, normally
    /// `i omega displacement`; no extra velocity or spin term is inferred.
    pub fn bearing_reactions(
        &self,
        displacement: &[C64],
        velocity: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<BearingReaction>, RotorError> {
        poll(&mut cancelled)?;
        if displacement.len() != 4
            || velocity.len() != 4
            || displacement
                .iter()
                .chain(velocity)
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(RotorError::Invalid(
                "four finite rotor coordinates and rates required",
            ));
        }
        let mut result = Vec::with_capacity(self.bearings.len());
        for (bearing, data) in self.bearings.iter().enumerate() {
            poll(&mut cancelled)?;
            let z = data.station_m;
            let u = [
                displacement[0] + displacement[3].scale(z),
                displacement[1] - displacement[2].scale(z),
            ];
            let v = [
                velocity[0] + velocity[3].scale(z),
                velocity[1] - velocity[2].scale(z),
            ];
            let force = std::array::from_fn(|i| {
                -(u[0].scale(data.stiffness_n_per_m[i][0])
                    + u[1].scale(data.stiffness_n_per_m[i][1])
                    + v[0].scale(data.damping_ns_per_m[i][0])
                    + v[1].scale(data.damping_ns_per_m[i][1]))
            });
            let moment = [force[1].scale(-z), force[0].scale(z)];
            let power = finite(0.5 * (force[0] * v[0].conj() + force[1] * v[1].conj()).re)?;
            if u.iter()
                .chain(&v)
                .chain(&force)
                .chain(&moment)
                .any(|x| !x.re.is_finite() || !x.im.is_finite())
            {
                return Err(RotorError::Invalid("unrepresentable bearing response"));
            }
            result.push(BearingReaction {
                bearing,
                displacement_m: u,
                force_n: force,
                moment_nm: moment,
                mean_power_on_rotor_w: power,
            });
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

fn point_map(z: f64) -> [[f64; 4]; 2] {
    [[1.0, 0.0, 0.0, z], [0.0, 1.0, -z, 0.0]]
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RotorError> {
    if cancelled() {
        Err(RotorError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, RotorError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RotorError::Invalid(
            "nonfinite rigid rotor coefficient or response",
        ))
    }
}
fn add_compensated(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), RotorError> {
    let next = finite(*sum + value)?;
    *correction = finite(
        *correction
            + if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}
