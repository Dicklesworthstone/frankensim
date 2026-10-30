//! Revolution-averaged rigid-blade loads in oblique flow.
//!
//! Each section sees `U_t = Ω r + V·t` and `U_a = V·n + v_i`, where
//! `V = V_hub + ω_tilt × r`. Native section thrust and drag resolve into
//! `dF = dT n - (dQ/r) t`; their actual arms supply the transverse hub
//! moment and axial shaft torque. This is the blade-element construction
//! described in NASA CR-182121, discussion of figures 93–95:
//! <https://ntrs.nasa.gov/citations/19880016976>.
//!
//! The model averages a full revolution, not a frozen blade phase. It keeps
//! the native radial quadrature and Estimated attached/separated section law.
//! Spanwise flow is reported but exerts no force in this two-dimensional
//! strip approximation. There is no flapping, compressibility, swirl,
//! skew-induced wake distribution, cyclic pitch or unsteady-airfoil memory.
//! Induction is supplied by the caller, separately from this blade-load law.

use super::dynamic_inflow::first_harmonic::HarmonicInduction;
use super::{MAX_BRAKING_INDUCTION, Refusal, Rotor, StationProblem, refuse};
use fs_math::det;

/// Hard work cap; construction requires 8–512 points in multiples of four.
pub const MAX_AZIMUTH_POINTS: usize = 512;

/// Physical spin direction relative to the independently specified thrust axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpinDirection {
    AlongThrust,
    AgainstThrust,
}

impl SpinDirection {
    fn sign(self) -> f64 {
        match self {
            Self::AlongThrust => 1.0,
            Self::AgainstThrust => -1.0,
        }
    }
}

/// Material disk frame. Rotate BOTH axis and reference in frame comparisons.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiskFrame {
    axis: [f64; 3],
    radial: [f64; 3],
    azimuth: [f64; 3],
    spin: SpinDirection,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn invalid(message: &str) -> Refusal {
    refuse(
        "oblique-rotor-invalid",
        message.into(),
        "supply a finite orthonormal disk frame, bounded quadrature and a normal-flow rotor state",
    )
}

impl DiskFrame {
    /// Admit orthonormal directions without silently changing their geometry.
    ///
    /// # Errors
    /// `oblique-rotor-invalid` for nonfinite, nonunit or nonorthogonal axes.
    pub fn new(
        thrust_axis: [f64; 3],
        radial_reference: [f64; 3],
        spin: SpinDirection,
    ) -> Result<Self, Refusal> {
        if !thrust_axis
            .iter()
            .chain(&radial_reference)
            .all(|x| x.is_finite())
            || (dot(thrust_axis, thrust_axis) - 1.0).abs() > 1e-12
            || (dot(radial_reference, radial_reference) - 1.0).abs() > 1e-12
            || dot(thrust_axis, radial_reference).abs() > 1e-12
        {
            return Err(invalid(
                "disk directions are not an orthonormal finite pair",
            ));
        }
        Ok(Self {
            axis: thrust_axis,
            radial: radial_reference,
            azimuth: cross(thrust_axis, radial_reference),
            spin,
        })
    }

    /// Positive thrust direction in the caller's frame.
    #[must_use]
    pub fn thrust_axis(self) -> [f64; 3] {
        self.axis
    }

    /// Nonspinning radial reference defining the cosine inflow mode.
    #[must_use]
    pub fn radial_reference(self) -> [f64; 3] {
        self.radial
    }

    /// Positive sine direction, independent of physical blade handedness.
    #[must_use]
    pub fn azimuth_reference(self) -> [f64; 3] {
        self.azimuth
    }

    /// Re-express the SAME induction field in another reference on this disk.
    /// A changed normal requires a moving-disk model rather than a gauge change.
    pub fn rebase_induction(
        self,
        other: Self,
        induction: HarmonicInduction,
    ) -> Result<HarmonicInduction, Refusal> {
        if (dot(self.axis, other.axis) - 1.0).abs() > 1e-12
            || !induction.modes().into_iter().all(f64::is_finite)
        {
            return Err(invalid("inflow rebase requires the same oriented disk"));
        }
        let gradient = core::array::from_fn(|j| {
            induction.first_mps[0] * self.radial[j] + induction.first_mps[1] * self.azimuth[j]
        });
        let result = HarmonicInduction {
            mean_mps: induction.mean_mps,
            first_mps: [dot(gradient, other.radial), dot(gradient, other.azimuth)],
        };
        if !result.modes().into_iter().all(f64::is_finite) {
            return Err(invalid("inflow rebase overflow"));
        }
        Ok(result)
    }

    /// Power-conjugate loads for `[mean, cosine, sine]` induction modes.
    /// The two moment loads are divided by outer disk radius and have units N.
    pub fn wake_force_n(self, loads: &ObliqueLoads, radius_m: f64) -> Result<[f64; 3], Refusal> {
        let force = [
            loads.thrust_n,
            -dot(loads.hub_moment_nm, self.azimuth) / radius_m,
            dot(loads.hub_moment_nm, self.radial) / radius_m,
        ];
        if !radius_m.is_finite() || radius_m <= 0.0 || !force.into_iter().all(f64::is_finite) {
            return Err(invalid("invalid radius or harmonic wake force"));
        }
        Ok(force)
    }
}

/// All loads are on the physical rotor, expressed in the input frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObliqueLoads {
    /// Full integrated aerodynamic force [N], including transverse force.
    pub force_n: [f64; 3],
    /// Transverse hub moment [N m]. EXCLUDES the axial shaft moment.
    pub hub_moment_nm: [f64; 3],
    /// Signed force along the positive thrust axis [N].
    pub thrust_n: f64,
    /// Required drive torque in the positive physical blade-speed convention.
    /// Actual aerodynamic shaft moment is `-spin_sign * torque_nm * axis`.
    pub torque_nm: f64,
    /// Positive profile-drag dissipation in the section-relative air [W].
    pub profile_dissipation_w: f64,
    /// Sum of section forces dotted with actual blade velocity relative to
    /// incident air, excluding the induced air velocity [W].
    pub kinematic_power_w: f64,
    /// Sum of section thrust times its actual induction [W], with its sign.
    /// In the uniform path this is exactly `thrust_n * induced_mps`.
    pub induced_work_w: f64,
    /// Hub force/tilt/absolute-spin power + induced work + profile loss [W].
    /// This is a section-law identity, not a whole-wake energy theorem.
    pub power_residual_w: f64,
    /// Largest ignored spanwise velocity [m/s], exposing the strip assumption.
    pub maximum_spanwise_speed_mps: f64,
    /// Lowest admitted tangential section speed [m/s], strictly positive.
    pub minimum_tangential_speed_mps: f64,
}

/// Prepared deterministic periodic quadrature. Evaluation allocates nothing;
/// storage and work scale with the admitted azimuth/station counts.
#[derive(Clone, Debug, PartialEq)]
pub struct AzimuthQuadrature {
    cos_sin: Vec<[f64; 2]>,
}

impl AzimuthQuadrature {
    /// Prepare a trapezoidal rule with explicit antipodal pairs.
    ///
    /// # Errors
    /// `oblique-rotor-invalid` outside 8–512 points or nonmultiples of four.
    pub fn new(points: usize) -> Result<Self, Refusal> {
        if !(8..=MAX_AZIMUTH_POINTS).contains(&points) || !points.is_multiple_of(4) {
            return Err(invalid(
                "azimuth count must be a multiple of four in [8, 512]",
            ));
        }
        let mut cos_sin = Vec::with_capacity(points);
        for k in 0..points / 4 {
            let angle = core::f64::consts::TAU * k as f64 / points as f64;
            let (c, s) = (det::cos(angle), det::sin(angle));
            cos_sin.extend_from_slice(&[[c, s], [-c, -s], [-s, c], [s, -c]]);
        }
        Ok(Self { cos_sin })
    }

    /// Fixed quadrature cost per radial station.
    #[must_use]
    pub fn points(&self) -> usize {
        self.cos_sin.len()
    }

    /// Integrate native blade sections over azimuth at the supplied uniform
    /// induction. `omega_rad_s` is positive ABSOLUTE spin about the physical
    /// spin axis. The carrier's axial rate is already included in that speed;
    /// only its transverse angular rate contributes additional section motion.
    /// Hub velocity is rotor material velocity relative to incident air.
    ///
    /// # Errors
    /// Native rotor/section admission; `oblique-rotor-invalid` for invalid
    /// inputs/overflow; `oblique-rotor-section-domain` for reversed local flow
    /// or angle outside the native unclipped section range `|alpha| <= 1.3`.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        rotor: &Rotor,
        density_kg_m3: f64,
        frame: DiskFrame,
        hub_velocity_mps: [f64; 3],
        carrier_angular_velocity_rad_s: [f64; 3],
        induced_mps: f64,
        omega_rad_s: f64,
    ) -> Result<ObliqueLoads, Refusal> {
        self.evaluate_modes(
            rotor,
            density_kg_m3,
            frame,
            hub_velocity_mps,
            carrier_angular_velocity_rad_s,
            HarmonicInduction {
                mean_mps: induced_mps,
                first_mps: [0.0; 2],
            },
            omega_rad_s,
            false,
        )
    }

    /// Evaluate first spatial harmonics at every real radial/azimuth point.
    /// Induced work is the actual sum `sum(dT * v_i)`, including the two
    /// transverse hub-moment ports. It is not thrust times mean induction.
    /// Exact zero harmonic amplitudes retain the uniform path's arithmetic.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_first_harmonic(
        &self,
        rotor: &Rotor,
        density_kg_m3: f64,
        frame: DiskFrame,
        hub_velocity_mps: [f64; 3],
        carrier_angular_velocity_rad_s: [f64; 3],
        induction: HarmonicInduction,
        omega_rad_s: f64,
    ) -> Result<ObliqueLoads, Refusal> {
        if !induction.modes().into_iter().all(f64::is_finite)
            || induction.first_mps[0].hypot(induction.first_mps[1])
                > dot(hub_velocity_mps, frame.axis) + induction.mean_mps
        {
            return Err(invalid(
                "invalid harmonic induction or reversed disk throughflow",
            ));
        }
        self.evaluate_modes(
            rotor,
            density_kg_m3,
            frame,
            hub_velocity_mps,
            carrier_angular_velocity_rad_s,
            induction,
            omega_rad_s,
            induction.first_mps != [0.0; 2],
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_modes(
        &self,
        rotor: &Rotor,
        density_kg_m3: f64,
        frame: DiskFrame,
        hub_velocity_mps: [f64; 3],
        carrier_angular_velocity_rad_s: [f64; 3],
        induction: HarmonicInduction,
        omega_rad_s: f64,
        distributed: bool,
    ) -> Result<ObliqueLoads, Refusal> {
        rotor.admit()?;
        let induced_mps = induction.mean_mps;
        let axial = dot(hub_velocity_mps, frame.axis);
        if !hub_velocity_mps
            .iter()
            .chain(&carrier_angular_velocity_rad_s)
            .chain(&[density_kg_m3, induced_mps, omega_rad_s, axial])
            .all(|x| x.is_finite())
            || density_kg_m3 <= 0.0
            || omega_rad_s <= 0.0
            || axial < 0.0
            || induced_mps < -MAX_BRAKING_INDUCTION * axial
        {
            return Err(invalid("invalid density, speed or normal-flow induction"));
        }
        let axial_rate = dot(carrier_angular_velocity_rad_s, frame.axis);
        let tilt = core::array::from_fn(|j| {
            carrier_angular_velocity_rad_s[j] - axial_rate * frame.axis[j]
        });
        let mut transverse_force = [0.0; 3];
        let mut hub_moment_nm = [0.0; 3];
        let mut thrust_n = 0.0;
        let mut torque_nm = 0.0;
        let mut profile_dissipation_w = 0.0;
        let mut kinematic_power_w = 0.0;
        let mut section_induced_work_w = 0.0;
        let mut maximum_spanwise_speed_mps: f64 = 0.0;
        let mut minimum_tangential_speed_mps = f64::INFINITY;
        let root_radius = rotor.stations[0].r_over_r * 0.8 * rotor.radius_m;
        for (index, station) in rotor.stations.iter().enumerate() {
            let radius = station.r_over_r * rotor.radius_m;
            let weight = rotor.station_width_m(index) / self.points() as f64;
            for &[c, s] in &self.cos_sin {
                let radial = core::array::from_fn(|j| c * frame.radial[j] + s * frame.azimuth[j]);
                let tangent = cross(frame.axis, radial).map(|x| frame.spin.sign() * x);
                let arm = radial.map(|x| radius * x);
                let tilt_velocity = cross(tilt, arm);
                let incident = core::array::from_fn(|j| hub_velocity_mps[j] + tilt_velocity[j]);
                let local_axial = dot(incident, frame.axis);
                let tangential = omega_rad_s * radius + dot(incident, tangent);
                let local_induced = if distributed {
                    induction.at(station.r_over_r * c, station.r_over_r * s)
                } else {
                    induced_mps
                };
                let throughflow = local_axial + local_induced;
                let alpha = station.beta_rad - det::atan2(throughflow, tangential);
                if ![tangential, throughflow, alpha]
                    .into_iter()
                    .all(f64::is_finite)
                {
                    return Err(invalid("section velocity overflow"));
                }
                if tangential <= 0.0 || throughflow < 0.0 || alpha.abs() > 1.3 {
                    return Err(refuse(
                        "oblique-rotor-section-domain",
                        format!(
                            "station {index} leaves positive tangential/normal flow or |alpha| <= 1.3"
                        ),
                        "use a smaller crossflow/tilt rate or an applicable reversed-flow section model",
                    ));
                }
                let loads = StationProblem {
                    rotor,
                    station,
                    index,
                    rho: density_kg_m3,
                    axial: local_axial,
                    tangential,
                    radius,
                    root_radius,
                }
                .blade_loads(local_induced)?;
                let d_thrust = loads.thrust * weight;
                let d_torque = loads.torque * weight;
                let d_tangential = d_torque / radius;
                let d_force =
                    core::array::from_fn(|j| d_thrust * frame.axis[j] - d_tangential * tangent[j]);
                let d_moment = cross(arm, frame.axis).map(|x| d_thrust * x);
                let material_velocity =
                    core::array::from_fn(|j| incident[j] + omega_rad_s * radius * tangent[j]);
                for j in 0..3 {
                    transverse_force[j] -= d_tangential * tangent[j];
                    hub_moment_nm[j] += d_moment[j];
                }
                thrust_n += d_thrust;
                if distributed {
                    section_induced_work_w += d_thrust * local_induced;
                }
                torque_nm += d_torque;
                let section_speed = (throughflow * throughflow + tangential * tangential).sqrt();
                profile_dissipation_w += loads.profile_drag * section_speed * weight;
                kinematic_power_w += dot(d_force, material_velocity);
                maximum_spanwise_speed_mps =
                    maximum_spanwise_speed_mps.max(dot(incident, radial).abs());
                minimum_tangential_speed_mps = minimum_tangential_speed_mps.min(tangential);
            }
        }
        let force_n = core::array::from_fn(|j| thrust_n * frame.axis[j] + transverse_force[j]);
        let induced_work_w = if distributed {
            section_induced_work_w
        } else {
            thrust_n * induced_mps
        };
        let power_residual_w = dot(force_n, hub_velocity_mps) + dot(hub_moment_nm, tilt)
            - torque_nm * omega_rad_s
            + induced_work_w
            + profile_dissipation_w;
        if !force_n
            .iter()
            .chain(&hub_moment_nm)
            .chain(&[
                thrust_n,
                torque_nm,
                profile_dissipation_w,
                kinematic_power_w,
                induced_work_w,
                power_residual_w,
                maximum_spanwise_speed_mps,
                minimum_tangential_speed_mps,
            ])
            .all(|x| x.is_finite())
        {
            return Err(invalid("integrated rotor force, moment or power overflow"));
        }
        Ok(ObliqueLoads {
            force_n,
            hub_moment_nm,
            thrust_n,
            torque_nm,
            profile_dissipation_w,
            kinematic_power_w,
            induced_work_w,
            power_residual_w,
            maximum_spanwise_speed_mps,
            minimum_tangential_speed_mps,
        })
    }
}
