//! Planar contact dynamics of the FS-1903 lever escapement.
//!
//! Three rigid bodies turn about fixed pivots in the tourbillon-cage frame:
//! the balance (angle θ, about `o_b`), the pallet lever (φ, about the origin) and
//! the escape wheel with the reflected inertia of the cage and train (ψ, about
//! `c_e`). They interact only through unilateral contacts between real outlines:
//!
//! - fifteen ratchet teeth (convex triangles) against the two pallet stones
//!   (parallelograms whose locking faces carry 12° of draw and whose impulse
//!   faces are inclined planes),
//! - the ruby impulse pin (a circle) against the three prongs of the fork slot,
//! - the lever against its two banking pins, and the guard dart against the
//!   safety roller outside the roller's passing crescent.
//!
//! Each step is a Moreau–Jean impulse step: smooth forces at the midpoint, a
//! simultaneous normal-impact LCP over all touching contacts solved by
//! `fs_contact::multicontact::solve_normal_impacts` (Delassus matrix
//! `W = J M⁻¹ Jᵀ`, Newton restitution), then Coulomb friction impulses from
//! `fs_contact::tangential::impact::solve_isotropic_coulomb_impulse` with a
//! Stribeck coefficient from `fs_tribo::FrictionLaw`. Locking, draw, unlocking,
//! impulse, drop, the beat, and therefore amplitude, rate and beat error are
//! outputs of this geometry, not inputs.
//!
//! While the impulse pin is far from the fork and the lever and wheel are at
//! rest, those two bodies sleep and the balance is advanced with the exact
//! damped-oscillator solution, so the expensive contact stepping runs only near
//! the line of centres (about an eighth of each beat).

use core::f64::consts::{PI, TAU};

use fs_contact::multicontact::{MultiImpactConfig, NormalContactImpact, solve_normal_impacts};
use fs_contact::tangential::impact::{CoulombImpulseInput, solve_isotropic_coulomb_impulse};
use fs_tribo::FrictionLaw;

const DEG: f64 = PI / 180.0;
const MM: f64 = 1.0e-3;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct V2 {
    pub x: f64,
    pub y: f64,
}

impl V2 {
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    #[must_use]
    pub fn polar(r: f64, a: f64) -> Self {
        Self::new(r * a.cos(), r * a.sin())
    }
    #[must_use]
    pub fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y)
    }
    #[must_use]
    pub fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y)
    }
    #[must_use]
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s)
    }
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y
    }
    #[must_use]
    pub fn cross(self, o: Self) -> f64 {
        self.x * o.y - self.y * o.x
    }
    #[must_use]
    pub fn len(self) -> f64 {
        // Not `hypot`: on wasm32 it is a soft-float libm call (through `fma`)
        // that dominated the whole engine's run time. Lengths here are
        // millimetre-scale, far from overflow or underflow.
        (self.x * self.x + self.y * self.y).sqrt()
    }
    #[must_use]
    pub fn rot(self, a: f64) -> Self {
        let (s, c) = a.sin_cos();
        Self::new(c * self.x - s * self.y, s * self.x + c * self.y)
    }
    #[must_use]
    pub fn perp(self) -> Self {
        Self::new(-self.y, self.x)
    }
    #[must_use]
    pub fn norm(self) -> Self {
        let l = self.len();
        if l > 0.0 { self.scale(1.0 / l) } else { Self::new(1.0, 0.0) }
    }
}

/// Make a polygon counter-clockwise (positive signed area).
fn ccw<const N: usize>(mut p: [V2; N]) -> [V2; N] {
    let mut area = 0.0;
    for i in 0..N {
        area += p[i].cross(p[(i + 1) % N]);
    }
    if area < 0.0 {
        p.reverse();
    }
    p
}

/// Closest point on segment `a`–`b` to `p`.
fn closest_on_segment(a: V2, b: V2, p: V2) -> V2 {
    let ab = b.sub(a);
    let t = (p.sub(a).dot(ab) / ab.dot(ab).max(1e-30)).clamp(0.0, 1.0);
    a.add(ab.scale(t))
}

/// Signed contact between convex CCW polygons `a` and `b`: gap (negative when
/// overlapping), unit normal pointing from `a` to `b`, and the witness points
/// on each body.
fn polygon_contact(a: &[V2], b: &[V2]) -> (f64, V2, V2, V2) {
    // Separating-axis test over both polygons' edge normals.
    let mut best_sep = f64::NEG_INFINITY;
    let mut best = (V2::default(), V2::default(), V2::default());
    for (poly, other, flip) in [(a, b, false), (b, a, true)] {
        let n = poly.len();
        for i in 0..n {
            let e0 = poly[i];
            let e1 = poly[(i + 1) % n];
            let out = e1.sub(e0).perp().scale(-1.0).norm(); // outward for CCW
            let (mut min_d, mut deepest) = (f64::INFINITY, other[0]);
            for &v in other {
                let d = out.dot(v.sub(e0));
                if d < min_d {
                    min_d = d;
                    deepest = v;
                }
            }
            if min_d > best_sep {
                best_sep = min_d;
                // normal from a to b
                let n_ab = if flip { out.scale(-1.0) } else { out };
                let on_face = deepest.sub(out.scale(min_d));
                best = if flip {
                    (n_ab, deepest, on_face)
                } else {
                    (n_ab, on_face, deepest)
                };
            }
        }
    }
    if best_sep <= 0.0 {
        return (best_sep, best.0, best.1, best.2);
    }
    // Separated: exact distance from vertex–segment pairs.
    let mut d_best = f64::INFINITY;
    let mut pa = V2::default();
    let mut pb = V2::default();
    for (poly, other, flip) in [(a, b, false), (b, a, true)] {
        let n = poly.len();
        for i in 0..n {
            let e0 = poly[i];
            let e1 = poly[(i + 1) % n];
            for &v in other {
                let c = closest_on_segment(e0, e1, v);
                let d = v.sub(c).len();
                if d < d_best {
                    d_best = d;
                    if flip {
                        pa = v;
                        pb = c;
                    } else {
                        pa = c;
                        pb = v;
                    }
                }
            }
        }
    }
    let n = pb.sub(pa).norm();
    (d_best, n, pa, pb)
}

/// Signed contact between convex CCW polygon `a` and a circle (`c`, `r`):
/// gap, normal from polygon to circle, witness points.
fn polygon_circle_contact(a: &[V2], c: V2, r: f64) -> (f64, V2, V2, V2) {
    let n = a.len();
    let mut s_max = f64::NEG_INFINITY;
    let mut n_max = V2::default();
    for i in 0..n {
        let e0 = a[i];
        let e1 = a[(i + 1) % n];
        let out = e1.sub(e0).perp().scale(-1.0).norm();
        let s = out.dot(c.sub(e0));
        if s > s_max {
            s_max = s;
            n_max = out;
        }
    }
    if s_max <= 0.0 {
        let pa = c.sub(n_max.scale(s_max));
        return (s_max - r, n_max, pa, c.sub(n_max.scale(r)));
    }
    let mut d_best = f64::INFINITY;
    let mut cp = a[0];
    for i in 0..n {
        let q = closest_on_segment(a[i], a[(i + 1) % n], c);
        let d = c.sub(q).len();
        if d < d_best {
            d_best = d;
            cp = q;
        }
    }
    let nrm = c.sub(cp).norm();
    (d_best - r, nrm, cp, c.sub(nrm.scale(r)))
}

/// Kinds of contact the escapement can make.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum ContactKind {
    EntryStoneTooth = 0,
    ExitStoneTooth = 1,
    PinFork = 2,
    Banking = 3,
    GuardRoller = 4,
}

/// One impact (a contact that closed with an approach speed above the
/// threshold) reported to the acoustics and the timegrapher.
#[derive(Clone, Copy, Debug)]
pub struct ContactImpact {
    pub kind: ContactKind,
    /// For stone–tooth contacts: true when the tooth landed on the locking face
    /// (a lock drop). For the pin: true on the unlocking (far) wall.
    pub is_lock: bool,
    pub time_s: f64,
    pub normal_impulse_n_s: f64,
    pub approach_speed_m_s: f64,
    pub dissipated_j: f64,
    /// World (cage-frame) contact point, metres.
    pub point: V2,
}

/// Fixed geometry and inertias of the escapement, designed from the classic
/// Swiss-lever proportions (15-tooth wheel, pallets embracing 2½ teeth,
/// 10.5° total lever angle, ~52° lift).
#[derive(Clone, Debug)]
pub struct EscapementDesign {
    pub n_teeth: usize,
    pub r_tip: f64,
    pub tooth_height: f64,
    pub tooth_lean: f64,
    pub tooth_back: f64,
    pub c_e: V2,
    pub o_b: V2,
    pub phi_bank: f64,
    pub phi_lock: f64,
    pub phi_run: f64,
    pub psi_drop: f64,
    pub draw: f64,
    pub half_lift: f64,
    pub entry_stone: [V2; 4],
    pub exit_stone: [V2; 4],
    pub r_roller: f64,
    pub r_pin: f64,
    pub fork: [[V2; 4]; 3],
    pub r_bank: f64,
    pub guard_gap_rad: f64,
    pub r_guard: f64,
    pub r_safety_roller: f64,
    pub crescent_half: f64,
    pub theta_zone: f64,
    pub i_lever: f64,
    pub i_escape_eff: f64,
    /// Tooth index offset so that ψ = 0 puts tooth 0 on the entry locking face.
    pub psi_tooth0: f64,
}

impl EscapementDesign {
    #[must_use]
    pub fn swiss_lever_pocket() -> Self {
        let n_teeth = 15;
        let pitch = TAU / n_teeth as f64;
        let r_tip = 2.40 * MM;
        let gamma = 1.25 * pitch; // half the angle embraced by the pallets (2½ teeth)
        let d_el = r_tip / gamma.cos();
        let c_e = V2::new(0.0, -d_el);
        let phi_bank = 5.25 * DEG;
        let phi_lock = 1.6 * DEG;
        let phi_run = 0.45 * DEG;
        let psi_drop = 1.8 * DEG;
        let draw = 12.0 * DEG;
        let psi_imp = 0.5 * pitch - psi_drop;
        let stone_len = 0.70 * MM;

        let entry_lock_polar = PI / 2.0 - gamma;
        let exit_lock_polar = PI / 2.0 + gamma;
        let make_stone = |lock_polar: f64, phi_unlock: f64, phi_letoff: f64| -> [V2; 4] {
            let l_world = c_e.add(V2::polar(r_tip, lock_polar));
            let p_world = c_e.add(V2::polar(r_tip, lock_polar + psi_imp));
            let l = l_world.rot(-phi_unlock);
            let p = p_world.rot(-phi_letoff);
            let f = V2::polar(1.0, lock_polar + draw).rot(-phi_unlock);
            ccw([l, p, p.add(f.scale(stone_len)), l.add(f.scale(stone_len))])
        };
        // Entry stone locks at −φ_b; exit at +φ_b (see module docs for the sense).
        let entry_stone = make_stone(
            entry_lock_polar,
            -phi_bank + phi_lock,
            phi_bank - phi_run,
        );
        let exit_stone = make_stone(
            exit_lock_polar,
            phi_bank - phi_lock,
            -phi_bank + phi_run,
        );

        // Fork and roller: the pin is aligned with the slot at θ = ±half_lift
        // when the lever rests on the corresponding banking.
        let half_lift = 26.0 * DEG;
        let r_roller = 0.75 * MM;
        let r_pin = 0.10 * MM;
        let d_bl =
            r_roller * half_lift.cos() + r_roller * half_lift.sin() / phi_bank.tan();
        let o_b = V2::new(0.0, d_bl);
        let slot_half = r_pin + 0.012 * MM;
        let prong = 0.22 * MM;
        let y_deep = d_bl - r_roller; // pin centre radius at θ = 0
        let y_bottom = y_deep - r_pin - 0.05 * MM;
        // Mouth below the pin path where it crosses the near prong (see design notes).
        let slope = (half_lift + phi_bank).tan(); // dρ/du of the pin path in the lever frame
        let y_align = (d_bl - r_roller * half_lift.cos()).hypot(r_roller * half_lift.sin());
        let y_mouth = y_align + slope * slot_half - r_pin * (1.0 + slope * slope).sqrt() - 0.01 * MM;
        let rect = |x0: f64, x1: f64, y0: f64, y1: f64| {
            ccw([V2::new(x0, y0), V2::new(x1, y0), V2::new(x1, y1), V2::new(x0, y1)])
        };
        let fork = [
            rect(-slot_half - prong, -slot_half, y_bottom, y_mouth),
            rect(slot_half, slot_half + prong, y_bottom, y_mouth),
            rect(-slot_half - prong, slot_half + prong, y_bottom - 0.15 * MM, y_bottom),
        ];

        // Inertias: steel lever ~10 mg; escape wheel plus 1/10² of the ~6e-9 kg m²
        // titanium cage that carries the escapement (the escape pinion rolls on the
        // fixed fourth wheel, ratio 10).
        let i_lever = 1.6e-11;
        let i_escape_eff = 2.2e-11 + 6.0e-9 / 100.0;

        Self {
            n_teeth,
            r_tip,
            tooth_height: 0.55 * MM,
            tooth_lean: 8.0 * DEG,
            tooth_back: 0.78 * pitch,
            c_e,
            o_b,
            phi_bank,
            phi_lock,
            phi_run,
            psi_drop,
            draw,
            half_lift,
            entry_stone,
            exit_stone,
            r_roller,
            r_pin,
            fork,
            r_bank: 1.2 * MM,
            guard_gap_rad: 0.30 * DEG,
            r_guard: 3.9 * MM,
            r_safety_roller: 0.45 * MM,
            crescent_half: half_lift + 8.0 * DEG,
            theta_zone: 42.0 * DEG,
            i_lever,
            i_escape_eff,
            psi_tooth0: entry_lock_polar,
        }
    }

    #[must_use]
    pub fn pitch(&self) -> f64 {
        TAU / self.n_teeth as f64
    }

    /// Tooth `i` outline in the cage frame for wheel angle `psi`.
    #[must_use]
    pub fn tooth(&self, i: usize, psi: f64) -> [V2; 3] {
        let a = psi + self.psi_tooth0 + i as f64 * self.pitch();
        let r_root = self.r_tip - self.tooth_height;
        ccw([
            self.c_e.add(V2::polar(self.r_tip, a)),
            self.c_e.add(V2::polar(r_root, a - self.tooth_lean)),
            self.c_e.add(V2::polar(r_root, a - self.tooth_back)),
        ])
    }

    /// Impulse-pin centre for balance angle `theta`.
    #[must_use]
    pub fn pin_center(&self, theta: f64) -> V2 {
        V2::new(self.r_roller * theta.sin(), self.o_b.y - self.r_roller * theta.cos())
    }
}

/// Smooth (non-contact) loads for one step, assembled by the horology layer
/// from the multiphysics state.
#[derive(Clone, Copy, Debug, Default)]
pub struct EscapementLoads {
    pub i_balance: f64,
    pub k_hairspring: f64,
    pub theta_eq: f64,
    /// Linear damping torque coefficient on the balance [N m s].
    pub balance_viscous: f64,
    /// Quadratic (form) drag coefficient [N m s²].
    pub balance_quadratic: f64,
    /// Coulomb friction torque on the balance pivots (or encasement) [N m].
    pub balance_coulomb: f64,
    /// Position-dependent external torque on the balance: amplitude and phase of
    /// a poise torque `a sin(θ + phase)` plus a constant.
    pub balance_poise_amp: f64,
    pub balance_poise_phase: f64,
    /// External torque on the lever (magnetic, out-of-poise) [N m].
    pub lever_torque: f64,
    pub lever_coulomb: f64,
    pub lever_viscous: f64,
    /// Drive torque at the escape wheel from the train [N m].
    pub escape_drive: f64,
    pub escape_viscous: f64,
    /// Coulomb brake on the escape wheel (pivots, frozen oil, solid gold) [N m].
    pub escape_coulomb: f64,
    /// Lubricated pallet friction (Stribeck) and contact restitutions.
    pub pallet_mu_static: f64,
    pub pallet_mu_kinetic: f64,
    pub pin_mu: f64,
}

#[derive(Clone, Copy, Debug)]
struct Contact {
    kind: ContactKind,
    gap: f64,
    j: [f64; 3],
    t: [f64; 3],
    restitution: f64,
    mu_static: f64,
    mu_kinetic: f64,
    is_lock: bool,
    point: V2,
    /// Identifier stable while the same feature pair stays in contact.
    key: u32,
}

/// Dynamic state of the escapement.
#[derive(Clone, Debug)]
pub struct Escapement {
    pub design: EscapementDesign,
    /// Generalised coordinates [θ balance, φ lever, ψ wheel] and velocities.
    pub q: [f64; 3],
    pub v: [f64; 3],
    pub time_s: f64,
    /// Lever and wheel are at rest against banking and locking face.
    pub asleep: bool,
    quiet_steps: u32,
    sleep_lever_torque: f64,
    sleep_drive: f64,
    /// Bitmask of the contact kinds that carried impulse in the last fine step.
    pub active_mask: u32,
    prev_keys: Vec<u32>,
    pub fine_step_s: f64,
    pub max_coarse_step_s: f64,
    /// Energy bookkeeping for efficiency: work done by the drive on the wheel and
    /// by the lever on the balance through the pin, since the last reset.
    pub drive_work_j: f64,
    pub impulse_work_on_balance_j: f64,
    pub solver_fallbacks: u64,
    contacts: Vec<Contact>,
    /// Balance inertia of the most recent loads (position projection weights).
    i_balance: f64,
    /// All three bodies at rest and held (stopped watch): integration is
    /// suspended until the loads change or the periodic probe runs.
    pub frozen: bool,
    still_steps: u32,
    frozen_signature: [f64; 8],
    frozen_probe_s: f64,
}

impl Escapement {
    #[must_use]
    pub fn new(design: EscapementDesign, theta0: f64, omega0: f64) -> Self {
        let mut e = Self {
            q: [theta0, -design.phi_bank, -0.6 * DEG],
            v: [omega0, 0.0, 0.0],
            design,
            time_s: 0.0,
            asleep: false,
            quiet_steps: 0,
            sleep_lever_torque: 0.0,
            sleep_drive: 0.0,
            active_mask: 0,
            prev_keys: Vec::with_capacity(8),
            fine_step_s: 1.5e-5,
            max_coarse_step_s: 8.0e-4,
            drive_work_j: 0.0,
            impulse_work_on_balance_j: 0.0,
            solver_fallbacks: 0,
            contacts: Vec::with_capacity(12),
            i_balance: 2.8e-9,
            frozen: false,
            still_steps: 0,
            frozen_signature: [0.0; 8],
            frozen_probe_s: 0.0,
        };
        e.project_positions();
        e
    }

    fn masses(loads: &EscapementLoads, d: &EscapementDesign) -> [f64; 3] {
        [loads.i_balance.max(1e-12), d.i_lever, d.i_escape_eff]
    }

    /// Collect every contact whose gap is below `margin` at configuration `q`.
    fn detect(&mut self, q: [f64; 3], margin: f64, loads: &EscapementLoads) {
        let d = &self.design;
        self.contacts.clear();
        let o_l = V2::default();
        let (theta, phi, psi) = (q[0], q[1], q[2]);
        let rot_lever = |p: V2| p.rot(phi);

        // Pallet stones against nearby teeth.
        let pitch = d.pitch();
        for (side, stone_l) in [(0u32, &d.entry_stone), (1u32, &d.exit_stone)] {
            let stone = stone_l.map(rot_lever);
            let center = stone.iter().fold(V2::default(), |a, &p| a.add(p)).scale(0.25);
            let stone_polar = center.sub(d.c_e);
            let a_stone = stone_polar.y.atan2(stone_polar.x);
            // Teeth whose tips lie within ±1.6 pitches of the stone (as seen from c_e).
            let base = a_stone - psi - d.psi_tooth0;
            let i_mid = (base / pitch).round() as i64;
            for di in -2..=1 {
                let i = (i_mid + di).rem_euclid(d.n_teeth as i64) as usize;
                let tooth = d.tooth(i, psi);
                let (gap, n, pa, pb) = polygon_contact(&stone, &tooth);
                if gap > margin {
                    continue;
                }
                let ra = pa.sub(o_l);
                let rb = pb.sub(d.c_e);
                let t = n.perp();
                // Lock if the touching stone face is the locking face (normal
                // mostly along the tooth's direction of travel).
                let travel = rb.perp().norm();
                let is_lock = n.dot(travel).abs() > 0.6;
                self.contacts.push(Contact {
                    kind: if side == 0 {
                        ContactKind::EntryStoneTooth
                    } else {
                        ContactKind::ExitStoneTooth
                    },
                    gap,
                    j: [0.0, -ra.cross(n), rb.cross(n)],
                    t: [0.0, -ra.cross(t), rb.cross(t)],
                    restitution: 0.10,
                    mu_static: loads.pallet_mu_static,
                    mu_kinetic: loads.pallet_mu_kinetic,
                    is_lock,
                    point: pb,
                    key: 100 + side * 32 + i as u32,
                });
            }
        }

        // Impulse pin against the fork prongs (only near the line of centres).
        let theta_rel = wrap_pi(theta);
        if theta_rel.abs() < d.crescent_half + 20.0 * DEG {
            let pc = d.pin_center(theta);
            for (k, prong) in d.fork.iter().enumerate() {
                let poly = prong.map(rot_lever);
                let (gap, n, pa, pb) = polygon_circle_contact(&poly, pc, d.r_pin);
                if gap > margin {
                    continue;
                }
                let ra = pa.sub(o_l);
                let rb = pb.sub(d.o_b);
                let t = n.perp();
                // The far prong (in the pin's direction of travel) is the unlocking face.
                let is_lock = n.dot(V2::new(theta_rel.signum(), 0.0).rot(phi)) > 0.0;
                self.contacts.push(Contact {
                    kind: ContactKind::PinFork,
                    gap,
                    j: [rb.cross(n), -ra.cross(n), 0.0],
                    t: [rb.cross(t), -ra.cross(t), 0.0],
                    restitution: 0.06,
                    mu_static: loads.pin_mu,
                    mu_kinetic: loads.pin_mu,
                    is_lock,
                    point: pb,
                    key: 200 + k as u32,
                });
            }
        }

        // Banking pins.
        for (sgn, key) in [(1.0, 300u32), (-1.0, 301u32)] {
            let gap = (d.phi_bank - sgn * phi) * d.r_bank;
            if gap <= margin {
                self.contacts.push(Contact {
                    kind: ContactKind::Banking,
                    gap,
                    j: [0.0, -sgn * d.r_bank, 0.0],
                    t: [0.0; 3],
                    restitution: 0.30,
                    mu_static: 0.0,
                    mu_kinetic: 0.0,
                    is_lock: false,
                    point: V2::polar(d.r_bank, PI / 2.0 + sgn * d.phi_bank).add(V2::new(0.0, -2.0 * d.r_bank)),
                    key,
                });
            }
        }

        // Guard dart against the safety roller, outside the passing crescent.
        if theta_rel.abs() > d.crescent_half {
            for (sgn, key) in [(1.0, 400u32), (-1.0, 401u32)] {
                // Lever on the +sgn side may not fall back past (φ_b − guard gap).
                if sgn * phi > 0.0 {
                    let gap = (sgn * phi - (d.phi_bank - d.guard_gap_rad)) * d.r_guard;
                    if gap <= margin {
                        self.contacts.push(Contact {
                            kind: ContactKind::GuardRoller,
                            gap,
                            j: [0.0, sgn * d.r_guard, 0.0],
                            t: [d.r_safety_roller, 0.0, 0.0],
                            restitution: 0.2,
                            mu_static: 0.15,
                            mu_kinetic: 0.15,
                            is_lock: false,
                            point: d.o_b.add(V2::new(0.0, -d.r_safety_roller)),
                            key,
                        });
                    }
                }
            }
        }
    }

    /// Conservative and driving generalised forces at configuration `q`. The
    /// dissipative terms are applied implicitly by [`Self::damp`] so that a
    /// balance in glassy oil or immersed in molten gold stays stable at the
    /// fine step.
    fn smooth_forces(q: [f64; 3], l: &EscapementLoads) -> [f64; 3] {
        let theta = q[0];
        let f_b = -l.k_hairspring * (theta - l.theta_eq)
            + l.balance_poise_amp * (theta + l.balance_poise_phase).sin();
        [f_b, l.lever_torque, l.escape_drive]
    }

    /// Backward-Euler viscous and (linearised) quadratic damping:
    /// `v ← v / (1 + h c_eff / m)`, unconditionally stable for any damping.
    fn damp(v: &mut [f64; 3], v_ref: [f64; 3], m: [f64; 3], l: &EscapementLoads, h: f64) {
        let c = [
            l.balance_viscous + l.balance_quadratic * v_ref[0].abs(),
            l.lever_viscous,
            l.escape_viscous,
        ];
        for k in 0..3 {
            v[k] /= 1.0 + h * c[k].max(0.0) / m[k];
        }
    }

    /// Apply implicit Coulomb friction (stick when the impulse capacity exceeds
    /// the momentum), returning the dissipated energy.
    fn coulomb_clamp(v: &mut f64, m: f64, torque: f64, h: f64) {
        if torque <= 0.0 {
            return;
        }
        let dv = torque * h / m;
        if v.abs() <= dv {
            *v = 0.0;
        } else {
            *v -= dv * v.signum();
        }
    }

    /// Advance by `dt`, choosing fine contact steps or coarse sleeping steps.
    pub fn advance(&mut self, dt: f64, loads: &EscapementLoads, impacts: &mut Vec<ContactImpact>) {
        let mut remaining = dt;
        let mut guard = 0u32;
        if self.frozen {
            let sig = Self::signature(loads);
            let changed = sig
                .iter()
                .zip(&self.frozen_signature)
                .any(|(a, b)| (a - b).abs() > 1.0e-3 * a.abs().max(b.abs()).max(1e-15));
            self.frozen_probe_s -= dt;
            if changed || self.frozen_probe_s <= 0.0 {
                self.frozen = false;
                self.still_steps = 0;
            } else {
                self.time_s += dt;
                return;
            }
        }
        while remaining > 1.0e-12 && guard < 2_000_000 {
            guard += 1;
            if self.frozen {
                self.time_s += remaining;
                break;
            }
            if self.asleep && self.should_wake(loads) {
                self.asleep = false;
                self.quiet_steps = 0;
            }
            if self.asleep {
                let h = self.coarse_step_size(remaining, loads);
                self.coarse_step(h, loads);
                remaining -= h;
            } else {
                let h = self.fine_step_s.min(remaining);
                self.fine_step(h, loads, impacts);
                remaining -= h;
                self.maybe_sleep(loads);
                self.maybe_freeze(loads);
            }
        }
    }

    fn signature(l: &EscapementLoads) -> [f64; 8] {
        [
            l.k_hairspring,
            l.theta_eq,
            l.balance_coulomb,
            l.balance_poise_amp,
            l.lever_torque,
            l.escape_drive,
            l.escape_coulomb,
            l.lever_coulomb,
        ]
    }

    /// Suspend integration once every body has been at rest for 6 ms of fine
    /// steps: a stopped watch then costs nothing until something changes. A
    /// probe re-runs the contact dynamics every 0.25 s.
    fn maybe_freeze(&mut self, loads: &EscapementLoads) {
        let still = self.v.iter().all(|w| w.abs() < 2.0e-4);
        if still {
            self.still_steps += 1;
        } else {
            self.still_steps = 0;
        }
        if self.still_steps > 400 {
            self.frozen = true;
            self.v = [0.0; 3];
            self.frozen_signature = Self::signature(loads);
            self.frozen_probe_s = 0.25;
            self.active_mask = 0;
        }
    }

    fn should_wake(&self, loads: &EscapementLoads) -> bool {
        let theta_rel = wrap_pi(self.q[0] - loads.theta_eq);
        if wrap_pi(self.q[0]).abs() <= self.design.theta_zone {
            return true;
        }
        let _ = theta_rel;
        let draw_scale = (loads.escape_drive.abs() * 0.15).max(1.0e-9);
        (loads.lever_torque - self.sleep_lever_torque).abs() > 0.2 * draw_scale
            || (loads.escape_drive - self.sleep_drive).abs() > 0.05 * self.sleep_drive.abs().max(1e-9)
    }

    fn maybe_sleep(&mut self, loads: &EscapementLoads) {
        let quiet = self.v[1].abs() < 2.0e-3 && self.v[2].abs() < 2.0e-3;
        let far = wrap_pi(self.q[0]).abs() > self.design.theta_zone + 2.0 * DEG;
        if quiet && far {
            self.quiet_steps += 1;
        } else {
            self.quiet_steps = 0;
        }
        if self.quiet_steps > 40 {
            self.asleep = true;
            self.v[1] = 0.0;
            self.v[2] = 0.0;
            self.sleep_lever_torque = loads.lever_torque;
            self.sleep_drive = loads.escape_drive;
            self.active_mask = 0;
            self.prev_keys.clear();
        }
    }

    fn coarse_step_size(&self, remaining: f64, _loads: &EscapementLoads) -> f64 {
        let theta = wrap_pi(self.q[0]);
        let w = self.v[0];
        let mut h = self.max_coarse_step_s.min(remaining);
        let approaching = theta * w < 0.0;
        if approaching && w.abs() > 1e-9 {
            let to_zone = (theta.abs() - self.design.theta_zone) / w.abs();
            if to_zone > 0.0 {
                h = h.min(to_zone + 1.0e-6);
            } else {
                h = h.min(self.fine_step_s);
            }
        }
        h.max(1.0e-7)
    }

    /// Balance-only step: exact damped linear oscillator for the hairspring and
    /// viscous terms, with the nonlinear torques applied as half kicks.
    fn coarse_step(&mut self, h: f64, l: &EscapementLoads) {
        let m = l.i_balance.max(1e-12);
        let kick = |v: &mut f64, q: f64, hk: f64| {
            let tq = l.balance_poise_amp * (q + l.balance_poise_phase).sin();
            *v += tq * hk / m;
            *v /= 1.0 + l.balance_quadratic * v.abs() * hk / m;
            Self::coulomb_clamp(v, m, l.balance_coulomb, hk);
        };
        let x0 = self.q[0] - l.theta_eq;
        // Balance held by static friction (pivots, frozen oil, solid gold).
        let static_torque = -l.k_hairspring * x0 + l.balance_poise_amp * (self.q[0] + l.balance_poise_phase).sin();
        if self.v[0].abs() < 1e-9 && static_torque.abs() <= l.balance_coulomb {
            self.v[0] = 0.0;
            self.time_s += h;
            return;
        }
        let mut v = self.v[0];
        kick(&mut v, self.q[0], 0.5 * h);
        let (x1, v1) = damped_oscillator(x0, v, l.k_hairspring / m, l.balance_viscous / m, h);
        self.q[0] = x1 + l.theta_eq;
        let mut v1m = v1;
        kick(&mut v1m, self.q[0], 0.5 * h);
        self.v[0] = v1m;
        self.time_s += h;
    }

    fn fine_step(&mut self, h: f64, l: &EscapementLoads, impacts: &mut Vec<ContactImpact>) {
        let m = Self::masses(l, &self.design);
        self.i_balance = m[0];
        let v0 = self.v;
        let q_mid = [
            self.q[0] + 0.5 * h * v0[0],
            self.q[1] + 0.5 * h * v0[1],
            self.q[2] + 0.5 * h * v0[2],
        ];
        let f = Self::smooth_forces(q_mid, l);
        let mut v = [
            v0[0] + h * f[0] / m[0],
            v0[1] + h * f[1] / m[1],
            v0[2] + h * f[2] / m[2],
        ];
        Self::damp(&mut v, v0, m, l, h);
        let v_free_balance = v[0];
        self.drive_work_j += l.escape_drive * 0.5 * (v0[2] + v[2]).max(0.0) * h;

        self.detect(q_mid, 2.0e-8, l);
        let mut mask = 0u32;
        let mut keys_now: Vec<u32> = Vec::with_capacity(self.contacts.len());
        if !self.contacts.is_empty() {
            let n = self.contacts.len();
            let jv = |c: &Contact, v: &[f64; 3]| c.j[0] * v[0] + c.j[1] * v[1] + c.j[2] * v[2];
            let mut delassus = vec![0.0; n * n];
            for (i, ci) in self.contacts.iter().enumerate() {
                for (k, ck) in self.contacts.iter().enumerate() {
                    delassus[i * n + k] =
                        ci.j[0] * ck.j[0] / m[0] + ci.j[1] * ck.j[1] / m[1] + ci.j[2] * ck.j[2] / m[2];
                }
            }
            let speeds_before: Vec<f64> = self.contacts.iter().map(|c| jv(c, &v)).collect();
            let impacts_in: Vec<NormalContactImpact> = self
                .contacts
                .iter()
                .zip(&speeds_before)
                .map(|(c, &s)| NormalContactImpact {
                    gap: 0.0,
                    activation_tolerance: 0.0,
                    // Push penetration out over a few steps (Baumgarte on the
                    // velocity target); restitution only for real impacts.
                    relative_normal_speed: s + 0.2 * c.gap.min(0.0) / h,
                    restitution: if s < -5.0e-3 { c.restitution } else { 0.0 },
                })
                .collect();
            let lambdas = match solve_normal_impacts(
                &impacts_in,
                &delassus,
                MultiImpactConfig::default(),
                || false,
            ) {
                Ok(r) => r.normal_impulses,
                Err(_) => {
                    self.solver_fallbacks += 1;
                    pgs_normal(&impacts_in, &delassus)
                }
            };
            for (c, &lam) in self.contacts.iter().zip(&lambdas) {
                if lam > 0.0 {
                    for k in 0..3 {
                        v[k] += c.j[k] * lam / m[k];
                    }
                }
            }
            // Coulomb friction impulses, capacity μ λ.
            for (c, &lam) in self.contacts.iter().zip(&lambdas) {
                if lam <= 0.0 || c.mu_static <= 0.0 {
                    continue;
                }
                let w_t = c.t[0] * c.t[0] / m[0] + c.t[1] * c.t[1] / m[1] + c.t[2] * c.t[2] / m[2];
                if w_t <= 0.0 {
                    continue;
                }
                let vt = c.t[0] * v[0] + c.t[1] * v[1] + c.t[2] * v[2];
                let mu = FrictionLaw::Stribeck {
                    static_mu: c.mu_static,
                    kinetic_mu: c.mu_kinetic,
                    characteristic_speed: 2.0e-3,
                    viscous_per_speed: 0.0,
                }
                .kinetic_coefficient(vt.abs())
                .unwrap_or(c.mu_kinetic);
                let law = FrictionLaw::Coulomb { static_mu: mu, kinetic_mu: mu };
                let res = solve_isotropic_coulomb_impulse(
                    CoulombImpulseInput {
                        normal_impulse: lam,
                        relative_tangent_velocity: [vt, 0.0],
                        inverse_tangent_mass: w_t,
                    },
                    &law,
                );
                let pt = match res {
                    Ok(r) => r.tangent_impulse[0],
                    Err(_) => (-vt / w_t).clamp(-mu * lam, mu * lam),
                };
                for k in 0..3 {
                    v[k] += c.t[k] * pt / m[k];
                }
            }
            // Report impacts and the active mask.
            for (i, c) in self.contacts.iter().enumerate() {
                let lam = lambdas[i];
                if lam <= 0.0 {
                    continue;
                }
                mask |= 1 << (c.kind as u32);
                if matches!(c.kind, ContactKind::ExitStoneTooth) {
                    mask |= 1 << 1;
                }
                keys_now.push(c.key);
                let s_before = speeds_before[i];
                let newly = !self.prev_keys.contains(&c.key);
                if newly && s_before < -2.0e-3 {
                    let s_after = jv(c, &v);
                    impacts.push(ContactImpact {
                        kind: c.kind,
                        is_lock: c.is_lock,
                        time_s: self.time_s + 0.5 * h,
                        normal_impulse_n_s: lam,
                        approach_speed_m_s: -s_before,
                        dissipated_j: (-0.5 * lam * (s_before + s_after)).max(0.0),
                        point: c.point,
                    });
                }
            }
        }
        // Work done on the balance by the fork through the pin (impulse delivery).
        if mask & (1 << ContactKind::PinFork as u32) != 0 {
            let dv_b = v[0] - v_free_balance;
            self.impulse_work_on_balance_j += m[0] * dv_b * 0.5 * (v[0] + v0[0]);
        }
        // Implicit Coulomb friction (pivots, frozen oil, solid gold) last.
        Self::coulomb_clamp(&mut v[0], m[0], l.balance_coulomb, h);
        Self::coulomb_clamp(&mut v[1], m[1], l.lever_coulomb, h);
        Self::coulomb_clamp(&mut v[2], m[2], l.escape_coulomb, h);

        for k in 0..3 {
            self.q[k] = q_mid[k] + 0.5 * h * v[k];
        }
        self.v = v;
        self.time_s += h;
        self.active_mask = mask;
        self.prev_keys = keys_now;
        // Free flight cannot leave more than a fraction of a micrometre of
        // overlap, which the next step's Baumgarte target removes; only steps
        // that carried contacts get the position projection.
        if !self.contacts.is_empty() {
            self.project_positions();
        }
    }

    /// Remove residual penetration by a few Gauss–Seidel position passes
    /// (non-energetic: velocities are untouched).
    fn project_positions(&mut self) {
        let loads = EscapementLoads {
            i_balance: self.i_balance,
            pallet_mu_static: 0.0,
            pallet_mu_kinetic: 0.0,
            pin_mu: 0.0,
            ..Default::default()
        };
        let m = [self.i_balance, self.design.i_lever, self.design.i_escape_eff];
        for _ in 0..3 {
            self.detect(self.q, -2.0e-7, &loads);
            if self.contacts.is_empty() {
                break;
            }
            for c in &self.contacts {
                let w = c.j[0] * c.j[0] / m[0] + c.j[1] * c.j[1] / m[1] + c.j[2] * c.j[2] / m[2];
                if w <= 0.0 {
                    continue;
                }
                let mu = -c.gap / w;
                for k in 0..3 {
                    self.q[k] += c.j[k] * mu / m[k];
                }
            }
        }
    }

    /// Balance angle wrapped to (−π, π] relative to the line of centres.
    #[must_use]
    pub fn balance_theta(&self) -> f64 {
        self.q[0]
    }
}

/// Exact solution of `x'' + c x' + k x = 0` over `h` from `(x0, v0)`, for
/// under-, critically and over-damped motion (`k`, `c` per unit inertia).
#[must_use]
pub fn damped_oscillator(x0: f64, v0: f64, k: f64, c: f64, h: f64) -> (f64, f64) {
    let k = k.max(0.0);
    let c = c.max(0.0);
    if k <= 0.0 {
        // Pure damping: v decays, x integrates it.
        if c * h < 1e-12 {
            return (x0 + v0 * h, v0);
        }
        let e = (-c * h).exp();
        return (x0 + v0 * (1.0 - e) / c, v0 * e);
    }
    let w0 = k.sqrt();
    let zeta = c / (2.0 * w0);
    let a = 0.5 * c;
    if zeta < 1.0 - 1e-9 {
        let wd = w0 * (1.0 - zeta * zeta).sqrt();
        let decay = (-a * h).exp();
        let (s, co) = (wd * h).sin_cos();
        let x1 = decay * (x0 * co + (v0 + a * x0) / wd * s);
        let v1 = decay * (v0 * co - (v0 * a + k * x0) / wd * s);
        (x1, v1)
    } else if zeta > 1.0 + 1e-9 {
        let disc = (a * a - k).sqrt();
        let r1 = -k / (a + disc); // slow root, free of cancellation
        let r2 = -a - disc; // fast root
        // x = A e^{r1 t} + B e^{r2 t}
        let b = (v0 - r1 * x0) / (r2 - r1);
        let aa = x0 - b;
        let e1 = (r1 * h).exp();
        let e2 = (r2 * h).exp();
        (aa * e1 + b * e2, aa * r1 * e1 + b * r2 * e2)
    } else {
        let e = (-w0 * h).exp();
        let b = v0 + w0 * x0;
        let x1 = (x0 + b * h) * e;
        let v1 = (b - w0 * (x0 + b * h)) * e;
        (x1, v1)
    }
}

fn wrap_pi(a: f64) -> f64 {
    let mut x = (a + PI) % TAU;
    if x < 0.0 {
        x += TAU;
    }
    x - PI
}

/// Projected Gauss–Seidel fallback for the normal LCP (only used if the
/// fs-contact active-set solver refuses a degenerate configuration).
fn pgs_normal(c: &[NormalContactImpact], w: &[f64]) -> Vec<f64> {
    let n = c.len();
    let mut lam = vec![0.0; n];
    let targets: Vec<f64> = c.iter().map(|x| -x.restitution * x.relative_normal_speed.min(0.0)).collect();
    for _ in 0..60 {
        for i in 0..n {
            let wii = w[i * n + i];
            if wii <= 0.0 {
                continue;
            }
            let mut s = c[i].relative_normal_speed;
            for k in 0..n {
                s += w[i * n + k] * lam[k];
            }
            lam[i] = (lam[i] + (targets[i] - s) / wii).max(0.0);
        }
    }
    lam
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn reference_loads() -> EscapementLoads {
        let i_b = 2.8e-9;
        let w = TAU * 2.5;
        EscapementLoads {
            i_balance: i_b,
            k_hairspring: i_b * w * w,
            theta_eq: 0.0,
            balance_viscous: i_b * w / 900.0,
            balance_quadratic: 0.0,
            balance_coulomb: 3.0e-9,
            balance_poise_amp: 0.0,
            balance_poise_phase: 0.0,
            lever_torque: 0.0,
            lever_coulomb: 0.0,
            lever_viscous: 1.0e-14,
            escape_drive: 1.4e-6,
            escape_viscous: 1.0e-14,
            escape_coulomb: 2.0e-8,
            pallet_mu_static: 0.18,
            pallet_mu_kinetic: 0.12,
            pin_mu: 0.15,
        }
    }

    #[test]
    fn design_locks_and_draws_at_rest() {
        let d = EscapementDesign::swiss_lever_pocket();
        let mut e = Escapement::new(d, 4.6, 0.0);
        let loads = reference_loads();
        let mut impacts = Vec::new();
        // Balance far from the fork: lever and wheel settle into lock.
        for _ in 0..2000 {
            e.fine_step(1.5e-5, &loads, &mut impacts);
        }
        assert!(e.v[2].abs() < 1e-2, "wheel should be locked, ω = {}", e.v[2]);
        assert!((e.q[1] + e.design.phi_bank).abs() < 0.3 * DEG, "lever drawn to banking: φ = {}", e.q[1] / DEG);
        assert!(e.active_mask & 1 != 0, "entry stone carries the lock: mask {:b}", e.active_mask);
    }

    #[test]
    fn free_running_escapement_beats_and_sustains_amplitude() {
        let d = EscapementDesign::swiss_lever_pocket();
        let mut e = Escapement::new(d, 270.0 * DEG, 0.0);
        let loads = reference_loads();
        let mut impacts = Vec::new();
        let mut peak = 0.0_f64;
        let mut peaks = Vec::new();
        let mut last_w = 0.0;
        for _ in 0..800 {
            e.advance(0.01, &loads, &mut impacts);
            peak = peak.max(e.q[0].abs());
            if last_w * e.v[0] < 0.0 {
                peaks.push(peak);
                peak = 0.0;
            }
            last_w = e.v[0];
        }
        let mut beats: Vec<f64> = Vec::new();
        for im in impacts.iter().filter(|i| {
            matches!(i.kind, ContactKind::EntryStoneTooth | ContactKind::ExitStoneTooth) && i.is_lock
        }) {
            if beats.last().is_none_or(|&t| im.time_s - t > 0.02) {
                beats.push(im.time_s);
            }
        }
        let n_locks = beats.len();
        // 8 s at 18 000 vph = 40 beats.
        assert!((36..=44).contains(&n_locks), "lock drops in 8 s: {n_locks}; fallbacks {}", e.solver_fallbacks);
        let amp = peaks.iter().rev().take(6).sum::<f64>() / 6.0 / DEG;
        assert!(amp > 200.0 && amp < 330.0, "amplitude {amp}°");
        // Wheel advanced half a tooth per beat.
        let turns = e.q[2] / (TAU / 30.0);
        assert!((turns - n_locks as f64).abs() <= 2.0, "ψ steps {turns} vs locks {n_locks}");
    }
    #[test]
    #[ignore = "diagnostic dump"]
    fn dump_geometry_and_first_steps() {
        let d = EscapementDesign::swiss_lever_pocket();
        let mm = |p: V2| format!("({:.4},{:.4})", p.x / MM, p.y / MM);
        println!("o_b {} c_e {}", mm(d.o_b), mm(d.c_e));
        println!("entry {:?}", d.entry_stone.map(|p| mm(p)));
        println!("exit  {:?}", d.exit_stone.map(|p| mm(p)));
        for (k, f) in d.fork.iter().enumerate() {
            println!("fork{k} {:?}", f.map(|p| mm(p)));
        }
        let mut e = Escapement::new(d.clone(), 4.6, 0.0);
        println!("after projection q = [{:.4}, {:.4}°, {:.4}°]", e.q[0], e.q[1] / DEG, e.q[2] / DEG);
        for i in 0..15 {
            let t = d.tooth(i, e.q[2]);
            let a = t[0].sub(d.c_e);
            let tip = t.iter().copied().fold(t[0], |b, p| if p.sub(d.c_e).len() > b.sub(d.c_e).len() { p } else { b });
            let at = tip.sub(d.c_e);
            if (at.y.atan2(at.x) / DEG - 90.0).abs() < 50.0 {
                println!("tooth {i} tip polar {:.2}° {}", at.y.atan2(at.x) / DEG, mm(tip));
            }
            let _ = a;
        }
        let loads = reference_loads();
        let mut impacts = Vec::new();
        for step in 0..400 {
            e.detect(e.q, 1.0e-6, &loads);
            if step % 20 == 0 || step < 6 {
                let cs: Vec<String> = e.contacts.iter().map(|c| format!("{:?}{} g={:.2}um J=[{:.2e},{:.2e},{:.2e}]", c.kind, if c.is_lock {"L"} else {""}, c.gap*1e6, c.j[0], c.j[1], c.j[2])).collect();
                println!("{step:4} φ={:.3}° ψ={:.3}° ω=[{:.3},{:.4},{:.4}] {:?}", e.q[1]/DEG, e.q[2]/DEG, e.v[0], e.v[1], e.v[2], cs);
            }
            e.fine_step(1.5e-5, &loads, &mut impacts);
        }
    }

    #[test]
    #[ignore = "diagnostic dump"]
    fn dump_free_run() {
        let d = EscapementDesign::swiss_lever_pocket();
        let mut e = Escapement::new(d, 270.0 * DEG, 0.0);
        let loads = reference_loads();
        let mut impacts = Vec::new();
        let mut peak = 0.0_f64;
        let mut last_w = 0.0;
        for k in 0..400 {
            let n0 = impacts.len();
            e.advance(0.01, &loads, &mut impacts);
            peak = peak.max(e.q[0].abs());
            if last_w * e.v[0] < 0.0 {
                println!("t={:.3} swing peak {:.1}° φ={:.2}° ψ={:.2}° asleep={} ", e.time_s, peak / DEG, e.q[1] / DEG, e.q[2] / DEG, e.asleep);
                peak = 0.0;
            }
            last_w = e.v[0];
            for im in &impacts[n0..] {
                if k < 120 {
                    println!("   t={:.4} {:?} lock={} J={:.2e} v={:.4} θ={:.1}° φ={:.2}°", im.time_s, im.kind, im.is_lock, im.normal_impulse_n_s, im.approach_speed_m_s, e.q[0] / DEG, e.q[1] / DEG);
                }
            }
        }
        let mut beats = Vec::new();
        for im in impacts.iter().filter(|i| matches!(i.kind, ContactKind::EntryStoneTooth | ContactKind::ExitStoneTooth) && i.is_lock) {
            if beats.last().map_or(true, |&t: &f64| im.time_s - t > 0.02) { beats.push(im.time_s); }
        }
        let iv: Vec<f64> = beats.windows(2).map(|w| w[1] - w[0]).collect();
        println!("beats {} intervals {:?}", beats.len(), &iv[iv.len().saturating_sub(8)..]);
        println!("psi {:.2}° fallbacks {}", e.q[2] / DEG, e.solver_fallbacks);
    }

}
