//! Complex continuation of the native uniaxial Maxwell eigensystem.
//! For q=(p,0,z), K=No²(I-aaᵀ)+Ne²aaᵀ and q×H=-K E.
//! Extraordinary dispersion is f z²+2gp z+h p²-No²Ne²=0. We retain
//! both stable quadratic roots, polish each by a certified double-double
//! Newton step and select the outgoing one by the certified sign of Im z,
//! then by the certified normal flux (`certified_root`); a principal complex
//! square root alone does not select it for a tilted lossy tensor. All fields
//! use the shared native eigenvectors,
//! row-scaled boundary solve and compensated physical-point propagator.

use super::super::super::kernel::{binary_scale, unit_axis};
use super::super::{
    Incidence, OrientedCrystalField, OrientedCrystalFieldSide,
    kernel::{
        self as transparent, FieldSum, WaveGeometry, cross, eigenvectors, incident_coordinates,
    },
    slab::kernel::{finite, flux, positive, solve_boundary},
};
use fs_math::{c64::C64, dd::Dd};

const C0: f64 = 299_792_458.0;
type Vector = [C64; 3];
type Result<T> = core::result::Result<T, &'static str>;

#[path = "kernel/absorption.rs"]
pub(super) mod absorption;

#[path = "kernel/box_integral.rs"]
mod box_integral;
pub(super) use box_integral::{
    admit_slab_box, coherent_box_absorption, coherent_box_requirements,
    coherent_slab_box_absorption, slab_box_requirements,
};

#[path = "kernel/certified_root.rs"]
mod certified_root;

#[path = "kernel/statistics.rs"]
mod statistics;
pub use statistics::{AbsorbingCrystalIlluminationPower, AbsorbingCrystalStatistics};
pub(super) use statistics::{
    IlluminationError, StatisticsSum, coherent_illumination_field, coherent_illumination_poynting,
    illumination_absorption, illumination_field_statistics, illumination_power,
    illumination_poynting, poll_illumination,
};

pub(super) fn slab_volume_m3(area: f64, start: f64, end: f64) -> Result<f64> {
    positive(area)?;
    if !start.is_finite() || !end.is_finite() || start < 0.0 || end <= start {
        return Err("crystal absorption needs a finite positive depth interval");
    }
    positive(area * (end - start))
}

/// One actual transmitted mode. An absorbing mode is an inhomogeneous wave;
/// its complex normal wavevector is not classified as a lossless bulk ray.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingCrystalMode {
    pub normal_wavenumber_rad_m: C64,
    pub electric_v_m: Vector,
    pub magnetic_a_m: Vector,
    /// This mode alone at the entrance. Cross-mode power is retained by the
    /// response and physical-point queries, not added into this number.
    pub normal_power_w_m2: f64,
}

/// Physical interface fields and power; transmitted power is the power
/// entering the crystal at z=0, before its subsequent volumetric absorption.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingCrystalResponse {
    /// Signed magnitude along `incidence_plane_axis()`, not necessarily kx.
    pub tangential_wavenumber_rad_m: f64,
    pub incident_normal_wavenumber_rad_m: f64,
    pub incident_electric_v_m: Vector,
    pub incident_magnetic_a_m: Vector,
    pub reflected_electric_v_m: Vector,
    pub reflected_magnetic_a_m: Vector,
    /// Ordinary then extraordinary; isotropic degeneracy uses TE then TM.
    pub transmitted: [AbsorbingCrystalMode; 2],
    pub incident_normal_power_w_m2: f64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
    omega: f64,
    axis: [f64; 3],
    loss_f_m: [f64; 2],
    incidence_plane_axis: [f64; 2],
}

/// Coherent peak-phasor E/H and energy transport, with material heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingCrystalField {
    pub field: OrientedCrystalField,
    /// omega/2 E*·[-Im(epsilon)]·E [W/m³]. This is zero in the lossless
    /// incident half-space and for pure evanescence in a lossless eigenaxis.
    pub absorbed_power_density_w_m3: f64,
}

impl AbsorbingCrystalResponse {
    pub(super) fn slab_absorbed_power_w(&self, area: f64, start: f64, end: f64) -> Result<f64> {
        absorption::integrated_heat(self, area, start, end)
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }

    /// Normalized transverse incidence direction in the original global
    /// interface frame. Multiply by the signed tangential wavenumber to
    /// obtain the physical [kx,ky] components.
    pub const fn incidence_plane_axis(&self) -> [f64; 2] {
        self.incidence_plane_axis
    }

    pub const fn principal_loss_permittivities_f_m(&self) -> [f64; 2] {
        self.loss_f_m
    }

    pub fn field(
        &self,
        side: OrientedCrystalFieldSide,
        point: [f64; 3],
    ) -> Result<AbsorbingCrystalField> {
        let field = self.electromagnetic_field(side, point)?;
        self.with_heating(side, field)
    }

    fn with_heating(
        &self,
        side: OrientedCrystalFieldSide,
        field: OrientedCrystalField,
    ) -> Result<AbsorbingCrystalField> {
        use OrientedCrystalFieldSide::*;
        let heat = if matches!(side, Ordinary | Extraordinary | Transmitted) {
            heating(field.electric_v_m, self.axis, self.loss_f_m, self.omega)?
        } else {
            0.0
        };
        Ok(AbsorbingCrystalField {
            field,
            absorbed_power_density_w_m3: heat,
        })
    }

    // Flux observations use this same coherent assembly without requiring
    // unrelated local dielectric work or electric squared moments to fit.
    fn electromagnetic_field(
        &self,
        side: OrientedCrystalFieldSide,
        point: [f64; 3],
    ) -> Result<OrientedCrystalField> {
        let mut sum = FieldSum::new(side, point)?;
        self.add_waves(side, &mut sum)?;
        sum.finish()
    }

    // All incident directions share one compensated, relative-phase sum.
    // Do not round each beam's propagated field before combining beams.
    fn add_waves(&self, side: OrientedCrystalFieldSide, sum: &mut FieldSum) -> Result<()> {
        use OrientedCrystalFieldSide::*;
        let mut add = |e: &Vector, h: &Vector, normal| {
            sum.add_wave(
                e,
                h,
                WaveGeometry {
                    tangential_wavenumber_rad_m: self.tangential_wavenumber_rad_m,
                    incidence_plane_axis: self.incidence_plane_axis,
                    normal_wavenumber_rad_m: normal,
                    anchor_z_m: 0.0,
                },
            )
        };
        if matches!(side, Incident | TotalIncidentSide) {
            add(
                &self.incident_electric_v_m,
                &self.incident_magnetic_a_m,
                C64::from_re(self.incident_normal_wavenumber_rad_m),
            )?;
        }
        if matches!(side, Reflected | TotalIncidentSide) {
            add(
                &self.reflected_electric_v_m,
                &self.reflected_magnetic_a_m,
                C64::from_re(-self.incident_normal_wavenumber_rad_m),
            )?;
        }
        for (m, mode) in self.transmitted.iter().enumerate() {
            if side == Transmitted || side == Ordinary && m == 0 || side == Extraordinary && m == 1
            {
                add(
                    &mode.electric_v_m,
                    &mode.magnetic_a_m,
                    mode.normal_wavenumber_rad_m,
                )?;
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
pub(super) struct Input {
    pub(super) omega: f64,
    pub(super) ni: f64,
    pub(super) indices: [C64; 2],
    pub(super) mui: f64,
    pub(super) mut_: f64,
    pub(super) axis: [f64; 3],
    pub(super) incidence: Incidence,
    pub(super) electric: [C64; 2],
    pub(super) loss_f_m: [f64; 2],
}

#[derive(Clone, Copy)]
pub(super) struct Mode {
    pub(super) z: C64,
    pub(super) e: Vector,
    pub(super) h: Vector,
}

/// Passive coordinate rotation around the unchanged interface normal. The
/// original global material axis and global Cartesian Jones input are both
/// expressed in the local x-z plane for the same native Maxwell solve.
pub(super) fn solve_in_plane(
    mut data: Input,
    plane_axis: [f64; 2],
) -> Result<AbsorbingCrystalResponse> {
    let p = unit_axis(plane_axis)
        .map_err(|_| "absorbing incidence-plane axis must be finite and nonzero")?;
    if p == [1.0, 0.0] {
        return solve(data);
    }
    let global_axis = transparent::unit(data.axis)?;
    let [x, y, z] = global_axis;
    data.axis = [x * p[0] + y * p[1], y * p[0] - x * p[1], z];
    let [ex, ey] = data.electric;
    data.electric = [
        ex.scale(p[0]) + ey.scale(p[1]),
        ey.scale(p[0]) - ex.scale(p[1]),
    ];
    let mut response = solve(data)?;
    // Same proper rotation as the transparent angular-spectrum owner: both
    // polar E and axial H transform alike. A sign/index pair is not clipped.
    let rotate = |v: Vector| -> Result<Vector> {
        let global = [
            v[0].scale(p[0]) - v[1].scale(p[1]),
            v[0].scale(p[1]) + v[1].scale(p[0]),
            v[2],
        ];
        if global.iter().any(|&v| !finite(v)) {
            return Err("absorbing crystal global field rotation is unrepresentable");
        }
        Ok(global)
    };
    response.incident_electric_v_m = rotate(response.incident_electric_v_m)?;
    response.incident_magnetic_a_m = rotate(response.incident_magnetic_a_m)?;
    response.reflected_electric_v_m = rotate(response.reflected_electric_v_m)?;
    response.reflected_magnetic_a_m = rotate(response.reflected_magnetic_a_m)?;
    for mode in &mut response.transmitted {
        mode.electric_v_m = rotate(mode.electric_v_m)?;
        mode.magnetic_a_m = rotate(mode.magnetic_a_m)?;
    }
    // Dielectric work and all propagated phasors use the original global
    // tensor frame, including the analytic ordinary/extraordinary cross term.
    response.axis = global_axis;
    response.incidence_plane_axis = p;
    Ok(response)
}

/// Unweighted outgoing Maxwell modes, shared by a half-space and a finite film.
#[derive(Clone, Copy)]
pub(super) struct Eigenbasis {
    pub(super) scale: f64,
    pub(super) i: f64,
    pub(super) p: f64,
    pub(super) qi: f64,
    pub(super) modes: [Mode; 2],
}

pub(super) fn eigenbasis(data: Input) -> Result<Eigenbasis> {
    let axis = transparent::unit(data.axis)?;
    positive(data.omega)?;
    positive(data.ni)?;
    positive(data.mui)?;
    positive(data.mut_)?;
    if data
        .indices
        .iter()
        .any(|&n| !finite(n) || n.re <= 0.0 || n.im > 0.0)
    {
        return Err("absorbing crystal requires passive principal indices");
    }
    if data.indices.iter().all(|n| n.im == 0.0) {
        let b = transparent::eigenbasis(transparent::Input {
            omega: data.omega,
            ni: data.ni,
            no: data.indices[0].re,
            ne: data.indices[1].re,
            mui: data.mui,
            mut_: data.mut_,
            axis: data.axis,
            incidence: data.incidence,
            electric: data.electric,
        })?;
        if b.modes
            .iter()
            .any(|m| m.kind == super::super::super::kernel::UniaxialWaveKind::Critical)
        {
            return Err("coalesced critical crystal modes require a generalized slab basis");
        }
        return Ok(Eigenbasis {
            scale: b.scale,
            i: b.i,
            p: b.p,
            qi: b.qi,
            modes: b.modes.map(|m| Mode {
                z: m.z,
                e: m.e,
                h: m.h,
            }),
        });
    }
    let scale = binary_scale(data.indices.iter().fold(data.ni, |s, n| s.max(n.abs())));
    let i = positive(data.ni / scale)?;
    let indices = data.indices.map(|n| C64::new(n.re / scale, n.im / scale));
    if indices
        .iter()
        .zip(data.indices)
        .any(|(n, original)| n.re == 0.0 || original.im != 0.0 && n.im == 0.0)
    {
        return Err("scaled crystal principal index or extinction is unrepresentable");
    }
    let [o2, e2] = indices.map(|n| C64::new((n.re + n.im) * (n.re - n.im), 2.0 * n.re * n.im));
    if [o2, e2]
        .iter()
        .zip(indices)
        .any(|(n2, n)| !finite(*n2) || *n2 == C64::ZERO || n.im != 0.0 && n2.im == 0.0)
    {
        return Err("scaled crystal constitutive tensor is unrepresentable");
    }
    let (p, p2, qi) = incident_coordinates(data.ni, scale, i, data.incidence)?;
    let od = Dd::from_f64(indices[0].re);
    let kd = Dd::from_f64(-indices[0].im);
    let ordinary_square = C64::new(((od - kd) * (od + kd) - p2).to_f64(), o2.im);
    let mut zo = ordinary_square.sqrt();
    if zo.im > 0.0 {
        zo = -zo;
    }
    let make = |z, ordinary| make_mode(z, p, axis, o2, e2, ordinary, indices[0] == indices[1]);
    let ordinary = make(zo, true)?;
    let extraordinary = if indices[0] == indices[1] {
        make(zo, false)?
    } else {
        let a2 = axis.map(|v| v * v);
        let norm = a2.iter().sum::<f64>();
        let f = (o2.scale(a2[0] + a2[1]) + e2.scale(a2[2])).scale(1.0 / norm);
        let g = (e2 - o2).scale(axis[0] * axis[2] / norm);
        let h = (o2.scale(a2[1] + a2[2]) + e2.scale(a2[0])).scale(1.0 / norm);
        let determinant = o2 * (o2.scale(a2[1]) + e2.scale(a2[0] + a2[2])).scale(1.0 / norm);
        let discriminant = f * o2 * e2 - determinant.scale(p2.to_f64());
        let root = discriminant.sqrt();
        let gp = g.scale(p);
        let plus = -gp + root;
        let minus = -gp - root;
        let numerator = if plus.abs() >= minus.abs() {
            plus
        } else {
            minus
        };
        if !finite(root)
            || root == C64::ZERO
            || !finite(f)
            || f == C64::ZERO
            || numerator == C64::ZERO
        {
            return Err("absorbing extraordinary dispersion is unresolved or coalesced");
        }
        let constant = h.scale(p2.to_f64()) - o2 * e2;
        let roots = [numerator / f, constant / numerator];
        // These double-precision roots only start the certified step: their Im carries about 1e-16 |z| of
        // rounding that the inputs do not.
        let dispersion = certified_root::Dispersion::new(indices, axis, p, p2);
        let mut candidates = [None; 2];
        let mut directions = [certified_root::Direction::Incoming; 2];
        for (slot, start) in roots.into_iter().enumerate() {
            let root = dispersion
                .refine(start)
                .ok_or(certified_root::NOT_CERTIFIED)?;
            candidates[slot] = Some(make(root.z, false)?);
            directions[slot] =
                certified_root::classify(root.z.im, root.im_bound(), || dispersion.flux(&root))?;
        }
        candidates[certified_root::select(directions)?]
            .ok_or("absorbing extraordinary mode has no resolved outgoing branch")?
    };
    Ok(Eigenbasis {
        scale,
        i,
        p,
        qi,
        modes: [ordinary, extraordinary],
    })
}

pub(super) fn solve(data: Input) -> Result<AbsorbingCrystalResponse> {
    let axis = transparent::unit(data.axis)?;
    positive(data.omega)?;
    positive(data.ni)?;
    positive(data.mui)?;
    positive(data.mut_)?;
    if data
        .indices
        .iter()
        .any(|&n| !finite(n) || n.re <= 0.0 || n.im > 0.0)
        || data.loss_f_m.iter().any(|v| !v.is_finite() || *v < 0.0)
        || data.electric.iter().any(|&v| !finite(v))
        || data
            .indices
            .iter()
            .zip(data.loss_f_m)
            .any(|(n, loss)| (n.im == 0.0) != (loss == 0.0))
    {
        return Err("absorbing crystal requires passive original axes and finite illumination");
    }
    if data.indices.iter().all(|n| n.im == 0.0) {
        // Preserve the established transparent arithmetic and exact critical
        // support when the original source supplies exactly zero extinction.
        let r = transparent::solve(transparent::Input {
            omega: data.omega,
            ni: data.ni,
            no: data.indices[0].re,
            ne: data.indices[1].re,
            mui: data.mui,
            mut_: data.mut_,
            axis: data.axis,
            incidence: data.incidence,
            electric: data.electric,
        })?;
        return Ok(AbsorbingCrystalResponse {
            tangential_wavenumber_rad_m: r.tangential_wavenumber_rad_m,
            incident_normal_wavenumber_rad_m: r.incident_normal_wavenumber_rad_m,
            incident_electric_v_m: r.incident_electric_v_m,
            incident_magnetic_a_m: r.incident_magnetic_a_m,
            reflected_electric_v_m: r.reflected_electric_v_m,
            reflected_magnetic_a_m: r.reflected_magnetic_a_m,
            transmitted: r.transmitted.map(|m| AbsorbingCrystalMode {
                normal_wavenumber_rad_m: m.normal_wavenumber_rad_m,
                electric_v_m: m.electric_v_m,
                magnetic_a_m: m.magnetic_a_m,
                normal_power_w_m2: m.normal_power_w_m2,
            }),
            incident_normal_power_w_m2: r.incident_normal_power_w_m2,
            reflected_power_fraction: r.reflected_power_fraction,
            transmitted_power_fraction: r.transmitted_power_fraction,
            omega: data.omega,
            axis,
            loss_f_m: data.loss_f_m,
            incidence_plane_axis: [1.0, 0.0],
        });
    }
    let Eigenbasis {
        scale,
        i,
        p,
        qi,
        mut modes,
    } = eigenbasis(data)?;
    let amplitude = positive(
        data.electric
            .iter()
            .fold(0.0_f64, |m, z| m.max(z.re.abs()).max(z.im.abs())),
    )?;
    let ei = data
        .electric
        .map(|z| C64::new(z.re / amplitude, z.im / amplitude));
    let ratio_mu = positive(data.mut_ / data.mui)?;
    let matrix = [
        modes.map(|m| m.h[1].scale(qi) + m.e[0].scale(ratio_mu * i * i)),
        modes.map(|m| -m.h[0] + m.e[1].scale(ratio_mu * qi)),
    ];
    let coefficients = solve_boundary(
        matrix,
        [
            ei[0].scale(2.0 * ratio_mu * i * i),
            ei[1].scale(2.0 * ratio_mu * qi),
        ],
    )?;
    for (mode, coefficient) in modes.iter_mut().zip(coefficients) {
        mode.e = mode.e.map(|v| v * coefficient);
        mode.h = mode.h.map(|v| v * coefficient);
    }
    let te = core::array::from_fn(|j| modes[0].e[j] + modes[1].e[j]);
    let th = core::array::from_fn(|j| modes[0].h[j] + modes[1].h[j]);
    let ie = [ei[0], ei[1], ei[0].scale(-p / qi)];
    let ih = cross([C64::from_re(p), C64::ZERO, C64::from_re(qi)], ie).map(|v| v.scale(ratio_mu));
    let re = [te[0] - ei[0], te[1] - ei[1], (te[0] - ei[0]).scale(p / qi)];
    let rh = cross([C64::from_re(p), C64::ZERO, C64::from_re(-qi)], re).map(|v| v.scale(ratio_mu));
    let hbound = ih
        .iter()
        .chain(&rh)
        .chain(&th)
        .fold(0.0_f64, |s, v| s.max(v.abs()));
    if (0..2).any(|j| (ih[j] + rh[j] - th[j]).abs() > 4096.0 * f64::EPSILON * hbound) {
        return Err("absorbing crystal tangential Maxwell continuity is unresolved");
    }
    let incoming = positive(flux(ie, ih))?;
    let reflected = -flux(re, rh) / incoming;
    let transmitted = flux(te, th) / incoming;
    if !reflected.is_finite()
        || !transmitted.is_finite()
        || transmitted < -1e-12
        || (reflected + transmitted - 1.0).abs() > 4096.0 * f64::EPSILON
    {
        return Err("absorbing crystal boundary normal-power balance failed");
    }
    let hscale = positive(scale / (data.mut_ * C0))?;
    let kscale = positive((data.omega / C0) * scale)?;
    let result = AbsorbingCrystalResponse {
        tangential_wavenumber_rad_m: kscale * p,
        incident_normal_wavenumber_rad_m: kscale * qi,
        incident_electric_v_m: ie.map(|v| v.scale(amplitude)),
        incident_magnetic_a_m: ih.map(|v| v.scale(amplitude * hscale)),
        reflected_electric_v_m: re.map(|v| v.scale(amplitude)),
        reflected_magnetic_a_m: rh.map(|v| v.scale(amplitude * hscale)),
        transmitted: modes.map(|m| AbsorbingCrystalMode {
            normal_wavenumber_rad_m: m.z.scale(kscale),
            electric_v_m: m.e.map(|v| v.scale(amplitude)),
            magnetic_a_m: m.h.map(|v| v.scale(amplitude * hscale)),
            normal_power_w_m2: flux(m.e, m.h) * amplitude * amplitude * hscale,
        }),
        incident_normal_power_w_m2: incoming * amplitude * amplitude * hscale,
        reflected_power_fraction: reflected,
        transmitted_power_fraction: transmitted.max(0.0),
        omega: data.omega,
        axis,
        loss_f_m: data.loss_f_m,
        incidence_plane_axis: [1.0, 0.0],
    };
    positive(result.incident_normal_power_w_m2)?;
    positive(result.incident_normal_wavenumber_rad_m)?;
    if !result.tangential_wavenumber_rad_m.is_finite()
        || result
            .incident_electric_v_m
            .iter()
            .chain(&result.incident_magnetic_a_m)
            .chain(&result.reflected_electric_v_m)
            .chain(&result.reflected_magnetic_a_m)
            .any(|&v| !finite(v))
        || result.transmitted.iter().zip(modes).any(|(m, original)| {
            !finite(m.normal_wavenumber_rad_m)
                || original.z.im != 0.0 && m.normal_wavenumber_rad_m.im == 0.0
                || !m.normal_power_w_m2.is_finite()
                || m.electric_v_m
                    .iter()
                    .chain(&m.magnetic_a_m)
                    .any(|&v| !finite(v))
        })
    {
        return Err("absorbing crystal output fields or decay are unrepresentable");
    }
    Ok(result)
}

fn make_mode(
    z: C64,
    p: f64,
    axis: [f64; 3],
    o2: C64,
    e2: C64,
    ordinary: bool,
    isotropic: bool,
) -> Result<Mode> {
    if !finite(z) {
        return Err("absorbing crystal normal wavevector is unrepresentable");
    }
    let norm = axis.iter().map(|v| v * v).sum::<f64>();
    let q = [C64::from_re(p), C64::ZERO, z];
    let (e, h) = eigenvectors(axis, q, ordinary, isotropic, |v| {
        let projection = (0..3).fold(C64::ZERO, |s, j| s + v[j].scale(axis[j] / norm));
        core::array::from_fn(|j| v[j] / o2 + projection * (e2.recip() - o2.recip()).scale(axis[j]))
    })?;
    let projection = (0..3).fold(C64::ZERO, |s, j| s + e[j].scale(axis[j] / norm));
    let displacement: Vector =
        core::array::from_fn(|j| e[j] * o2 + projection * (e2 - o2).scale(axis[j]));
    let faraday = cross(q, e);
    let ampere = cross(q, h);
    let bound = faraday
        .iter()
        .chain(&ampere)
        .chain(&h)
        .chain(&displacement)
        .fold(0.0_f64, |s, v| s.max(v.abs()));
    if !bound.is_finite()
        || bound <= 0.0
        || (0..3).any(|j| {
            (faraday[j] - h[j]).abs() > 8192.0 * f64::EPSILON * bound
                || (ampere[j] + displacement[j]).abs() > 8192.0 * f64::EPSILON * bound
        })
    {
        return Err("absorbing crystal eigenmode Maxwell residual failed");
    }
    Ok(Mode { z, e, h })
}

pub(super) fn heating(electric: Vector, axis: [f64; 3], loss: [f64; 2], omega: f64) -> Result<f64> {
    let a = axis.map(Dd::from_f64);
    let norm = a.iter().fold(Dd::ZERO, |s, &v| s + v * v);
    let dot = [false, true].map(|imaginary| {
        (0..3).fold(Dd::ZERO, |s, j| {
            s + a[j]
                * Dd::from_f64(if imaginary {
                    electric[j].im
                } else {
                    electric[j].re
                })
        })
    });
    let parallel = (dot[0] * dot[0] + dot[1] * dot[1]) / norm;
    let perpendicular = (0..3).fold(Dd::ZERO, |s, j| {
        let re = Dd::from_f64(electric[j].re) - a[j] * dot[0] / norm;
        let im = Dd::from_f64(electric[j].im) - a[j] * dot[1] / norm;
        s + re * re + im * im
    });
    // Positive weighted squares avoid subtracting almost equal total and
    // parallel intensities in a strongly dichroic material.
    let heat = ((perpendicular * Dd::from_f64(loss[0]) + parallel * Dd::from_f64(loss[1]))
        * Dd::from_f64(0.5 * omega))
    .to_f64();
    if !heat.is_finite() || heat < 0.0 {
        return Err("anisotropic absorbed-power density is unrepresentable");
    }
    Ok(heat)
}

#[cfg(test)]
#[path = "kernel/tests.rs"]
mod tests;
