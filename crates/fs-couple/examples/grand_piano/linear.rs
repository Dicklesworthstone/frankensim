//! Reciprocal moving-boundary string/soundboard mechanics.
//!
//! With y(x)=sum(phi_n(x) q_n)+x/L*b, the string has fixed-interface
//! modal mass 1, cross inertia beta_n=sqrt(2 mu L)(-1)^(n+1)/(n pi),
//! and endpoint mass mu L/3. Completing the kinetic square with
//! z_n=q_n+beta_n*b leaves positive residual endpoint mass
//! mu L/3-sum(beta_n^2). fs-modal mass-normalizes the small loaded board.
//!
//! The resulting potential is sum(w_n^2 (z_n-beta_n*b)^2)/2 plus the
//! bare board and T/L endpoint stiffness. Thus string and board forces
//! are RECIPROCAL, not a gyroscopic/velocity-coupling substitute.
//!
//! Existing fs-couple exact-ZOH transitions advance each diagonal oscillator.
//! The cross-potential uses average-displacement forces; a small Schur solve
//! closes their work exactly. Symmetric velocity flows apply physical wood
//! and moving-boundary string damping, including reciprocal bridge reactions.
//! The coupling is second-order consistent, NOT the exact full propagator.
//! No new oscillator, eigensolver or matrix factorization is implemented here.
//!
//! All unison members and duplex segments of a course/polarization share its
//! bridge shape. Optional orthogonal transverse motion uses the SAME bank,
//! board coordinates and energy, not a second independent piano.
//! Reduce partial forces BEFORE projecting to the board, and project board
//! displacement ONCE per course. This exact reassociation changes roundoff,
//! not the retained model: O(partials + courses*board_modes + board_modes^2)
//! stepping replaces repeated O(partials*board_modes) projections. Silent
//! strings still participate; no voice stealing or sympathetic-tail cutoff.

use std::f64::consts::{PI, TAU};
use fs_couple::modal_acoustic_time::{ModalAcousticMode,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::{c64::C64, det};
use super::geometry::Course;
use super::steinway_scale;

/// Spatial damper pads in these same mass-normalized moving-bridge coordinates.
#[path = "dampers.rs"]
pub mod dampers;
/// Independent physical contact sites on a supplied longitudinal hammer face.
#[path = "hammer_footprint.rs"]
pub mod hammer_footprint;
/// Geometric string extension in the SAME moving-bridge coordinates.
#[path = "string_stretching.rs"]
pub mod string_stretching;
/// Physical velocity and strain-rate loss with reciprocal bridge reactions.
#[path = "string_damping.rs"]
pub mod string_damping;

/// Prepared reduced-board capacity, shared with geometry/CSV/audio consumers.
/// This is a memory/work ceiling, not a claim that every size is real-time.
pub const MAX_BOARD_MODES: usize = 128;
/// A bass string may need hundreds of physical partials before reaching the
/// output band. The frequency ceiling still stops retention; this is a per-
/// string ceiling, not an instruction to allocate inaudible oscillators.
pub use super::geometry::MAX_STRING_MODES;

#[derive(Clone, Debug)]
pub struct BoardMode {
    pub frequency_hz: f64,
    pub damping_ratio: f64,
    /// Mass-normalized displacement at each key's bridge, 1/sqrt(kg).
    pub bridge: [f64; 88],
    /// Integral of the surface mode, m^2/sqrt(kg); not a pressure transfer.
    pub volume: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct StringMode {
    pub omega: f64,
    pub beta: f64,
    pub a: f64,
    pub hammer_shape: f64,
    pub damper_shape: f64,
    pub string: usize,
}

#[derive(Clone, Debug)]
pub struct StringPort {
    pub course: usize,
    /// Physical unison member, independent of contact-site subdivision.
    pub member: usize,
    /// 0: original hammer/vertical plane; 1: orthogonal transverse plane.
    pub polarization: usize,
    /// Segment identity must not be inferred from the absence of hammer contact.
    pub duplex: bool,
    pub modes: std::ops::Range<usize>,
    pub bridge: Vec<f64>,
    pub hammer_lift: f64,
    pub damper_lift: f64,
    pub contact: Option<usize>,
}

#[derive(Clone, Copy, Debug, Default)]
struct Transition { qq: f64, qv: f64, vq: f64, vv: f64, bq: f64, bv: f64 }

/// Cold lowering of the EXISTING public exact-ZOH implementation to a reusable
/// linear map. Unit-energy basis probes avoid extreme-frequency energy limits.
fn transition_chunk(rate: u32, modes: &[ModalAcousticMode]) -> Result<Vec<Transition>, String> {
    let model = ModalAcousticTimeModel::try_new(rate, modes.to_vec(),
        ModalAcousticTimeBudget::audible_reference()).map_err(|e| e.to_string())?;
    let mut result = Vec::with_capacity(modes.len());
    for i in 0..modes.len() {
        let (qq, qv, vq, vv, bq, bv) = model.probe_linear_map(i);
        if !bq.is_finite() || bq <= 0.0 { return Err("nonpositive ZOH displacement compliance".into()); }
        result.push(Transition { qq, qv, vq, vv, bq, bv });
    }
    Ok(result)
}

fn transitions(rate: u32, modes: &[ModalAcousticMode]) -> Result<Vec<Transition>, String> {
    let mut result = Vec::with_capacity(modes.len());
    for chunk in modes.chunks(4096) { result.extend(transition_chunk(rate, chunk)?); }
    Ok(result)
}

fn identity(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n { a[i * n + i] = 1.0; }
    a
}

/// All dense inversion is cold, through the existing symmetric eigenfacility.
fn inverse_spd(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let modes = fs_modal::eigh_gen_dense(a, &identity(n), n).map_err(|e| e.to_string())?;
    let mut inverse = vec![0.0; n * n];
    for m in modes {
        if !m.lambda.is_finite() || m.lambda <= 0.0 || m.residual > 1e-7 * m.lambda {
            return Err("bridge Schur complement is not resolved positive definite".into());
        }
        for i in 0..n { for j in 0..n { inverse[i*n+j] += m.phi[i]*m.phi[j]/m.lambda; } }
    }
    Ok(inverse)
}

/// Exact half-step velocity flow of the PHYSICAL bare-board damping after
/// string mass loading: C_loaded = Phi^T diag(2*zeta_i*omega_i) Phi.
/// The original per-mode loss is retained, including its off-diagonal terms.
/// Eigenanalysis and exponentials are cold, using the existing fs-modal owner.
/// Splitting this dissipative flow around the conservative bridge step is
/// second-order, not the exact full damped coupled propagator.
fn board_damping_flow(board: &[BoardMode], basis: &[f64], rate: u32)
    -> Result<Vec<f64>, String> {
    let r=board.len();
    let mut c=vec![0.0;r*r];
    for i in 0..r { for j in i..r {
        let value: f64=board.iter().enumerate().map(|(a,b)|
            2.0*b.damping_ratio*TAU*b.frequency_hz*basis[a*r+i]*basis[a*r+j]).sum();
        c[i*r+j]=value;c[j*r+i]=value;
    } }
    board_damping_matrix_flow(&c,r,rate)
}

fn board_damping_spectrum(c: &[f64], r: usize) -> Result<Vec<fs_modal::ModePair>, String> {
    if c.iter().any(|x|!x.is_finite()) {return Err("loaded board damping overflow".into());}
    let scale=c.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
    if !scale.is_finite() {return Err("loaded board damping overflow".into());}
    if scale==0.0 {return Ok(Vec::new());}
    let eigen=fs_modal::eigh_gen_dense(c,&identity(r),r).map_err(|e|e.to_string())?;
    let tolerance=1e-10*scale*r as f64;
    for pair in &eigen {
        if !pair.lambda.is_finite() || pair.lambda < -tolerance
            || !pair.residual.is_finite() || pair.residual > tolerance {
            return Err("physical board damping is not resolved positive semidefinite".into());
        }
    }
    Ok(eigen)
}

fn board_damping_matrix_flow(c: &[f64], r: usize, rate: u32) -> Result<Vec<f64>, String> {
    let eigen=board_damping_spectrum(c,r)?;
    if eigen.is_empty() {return Ok(identity(r));}
    let mut flow=vec![0.0;r*r];
    for pair in eigen {
        // A zero-loss direction can have a roundoff-negative eigenvalue.
        // Only the scale-aware numerical nullspace is projected to zero.
        let decay=det::exp(-0.5*pair.lambda.max(0.0)/f64::from(rate));
        for i in 0..r { for j in 0..r {flow[i*r+j]+=decay*pair.phi[i]*pair.phi[j];} }
    }
    Ok(flow)
}

/// No allocation. Workless damping changes kinetic energy only. Reporting
/// the actual norm difference closes the same floating-point energy ledger.
fn damp_board(flow: &[f64], from: &[f64], to: &mut [f64]) -> f64 {
    let r=from.len();
    for (i,out) in to.iter_mut().enumerate() {
        *out=flow[i*r..(i+1)*r].iter().zip(from).map(|(a,v)|a*v).sum();
    }
    0.5*from.iter().zip(to).map(|(a,b)|a*a-*b * *b).sum::<f64>()
}

/// `Bank::gpu_image`: the operators of a contact-free step. `transition[k]` is [qq, qv, vq, vv, bq, bv]
/// of coordinate k (string modes, then board coordinates); `groups` are the string ranges sharing one
/// bridge row (the first string's); `omega2` is the diagonal stiffness used in the modal-loss ledger.
#[derive(Clone, Debug)]
#[allow(dead_code)] // used by apps/instrument-lab/render-gpu-piano
pub struct BankImage {
    pub n: usize,
    pub r: usize,
    pub transition: Vec<[f64; 6]>,
    pub modes: Vec<StringMode>,
    pub omega2: Vec<f64>,
    pub strings: Vec<StringPort>,
    pub groups: Vec<std::ops::Range<usize>>,
    pub schur_inverse: Vec<f64>,
    pub board_half_damping: Option<Vec<f64>>,
    pub string_loss: Option<string_damping::LossImage>,
    pub contact_strings: Vec<usize>,
    pub q: Vec<f64>,
    pub v: Vec<f64>,
}

pub struct Bank {
    pub strings: Vec<StringPort>,
    pub modes: Vec<StringMode>,
    pub contact_strings: Vec<usize>,
    pub board_count: usize,
    pub q: Vec<f64>,
    pub v: Vec<f64>,
    pub next_q: Vec<f64>,
    pub next_v: Vec<f64>,
    /// Symmetric contact displacement compliance, m/N, row-major.
    pub contact_compliance: Vec<f64>,
    pub free_contact: Vec<f64>,
    pub rate: u32,
    footprints: Option<hammer_footprint::Prepared>,
    stretching: Option<string_stretching::Prepared>,
    /// Held work-conjugate loads in the original mass-normalized coordinates.
    /// Prepared contact owners update this scratch inside the same trial solve.
    external_force: Option<Vec<f64>>,
    pub omitted_duplex_modes: usize,
    /// Segment/direction count whose partial budget omits in-band vibration.
    pub capped_string_segments: usize,
    /// Contiguous strings sharing one physical course bridge; prepared cold.
    groups: Vec<std::ops::Range<usize>>,
    transition: Vec<Transition>,
    diagonal_omega2: Vec<f64>,
    pub last_modal_loss_j: f64,
    physical_board_k: Vec<f64>,
    board_half_damping: Option<Vec<f64>>,
    board_damping_enabled: bool,
    bare_board_damping_configured: bool,
    pub(super) string_loss: string_damping::Prepared,
    damped_v: Vec<f64>,
    board_pre_loss: f64,
    board_volume: Vec<f64>,
    /// Columns map loaded coordinates to the supplied bare-board coordinates.
    board_basis: Vec<f64>,
    schur_inverse: Vec<f64>,
    /// Row c is the prepared S^-1 W^T column for contact c.
    contact_response: Vec<f64>,
    #[cfg(test)]
    contact_board: Vec<f64>,
    free_q: Vec<f64>,
    free_v: Vec<f64>,
    r_string: Vec<f64>,
    board_rhs: Vec<f64>,
    board_end: Vec<f64>,
}

impl Bank {
    /// `rate` is the mechanics rate; `band_hz` is the OUTPUT band, so
    /// oversampling never silently admits inaudible/aliased retained modes.
    pub fn new(courses: &[Course], board: &[BoardMode], rate: u32, band_hz: f64,
        max_modes: usize, damping: bool) -> Result<Self, String> {
        Self::new_with_transverse_bridge(courses,board,rate,band_hz,max_modes,damping,None)
    }

    /// Optional orthogonal transverse bridge rows, course-major, in the SAME
    /// bare-board basis as BoardMode::bridge. The geometry consumer must project
    /// a physical unit direction orthogonal to the hammer plane and string axis.
    /// Both polarizations retain every speaking/unison/duplex mode admitted by
    /// the original per-span frequency budget, with the same physical T, EI, mu.
    /// Only the original vertical plane is directly struck. There is ONE hammer
    /// per course, not twice the contact area, force or launch energy. Each
    /// direction contributes its own endpoint kinetic/potential quadratic form;
    /// this is vector motion of one string mass, not twice its scalar mass.
    /// None preserves the original arithmetic and retained coordinates.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_transverse_bridge(courses: &[Course], board: &[BoardMode], rate: u32,
        band_hz: f64, max_modes: usize, damping: bool, secondary: Option<&[Vec<f64>]>)
        -> Result<Self, String> {
        Self::new_with_string_damping(courses, board, rate, band_hz, max_modes,
            damping, secondary, false)
    }

    /// Use RT-0425's per-note intrinsic string loss on the existing scalar
    /// stiff-string modes. The source's transverse rate and strain-rate terms
    /// give fixed-bridge zeta = R_u/omega + eta_u*omega. With a moving bridge,
    /// apply their physical velocity quadratic form, including chord traction.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_string_damping(courses: &[Course], board: &[BoardMode], rate: u32,
        band_hz: f64, max_modes: usize, damping: bool, secondary: Option<&[Vec<f64>]>,
        rt0425_string_damping: bool) -> Result<Self, String> {
        let r = board.len();
        if courses.is_empty() || courses.len() > 88 || !(1..=MAX_BOARD_MODES).contains(&r)
            || !(1..=MAX_STRING_MODES).contains(&max_modes) || rate < 8_000
            || !band_hz.is_finite() || band_hz <= 0.0 || band_hz > 0.45*f64::from(rate) {
            return Err("invalid course, modal, frequency or sample-rate budget".into());
        }
        if secondary.is_some_and(|rows| rows.len()!=courses.len()
            || rows.iter().any(|row|row.len()!=r || row.iter().any(|x|!x.is_finite()))) {
            return Err("transverse bridge rows must cover every course and bare-board mode".into());
        }
        let planes=if secondary.is_some() {2}else{1};
        for b in board {
            if !b.frequency_hz.is_finite() || b.frequency_hz <= 0.0
                || !b.damping_ratio.is_finite() || b.damping_ratio < 0.0
                || !b.volume.is_finite() || b.bridge.iter().any(|x| !x.is_finite()) {
                return Err("invalid mass-normalized soundboard mode".into());
            }
        }
        let mut strings = Vec::new();
        let mut groups = Vec::with_capacity(courses.len()*planes);
        let mut modes = Vec::new();
        let mut contact_strings = Vec::new();
        let mut oscillator = Vec::new();
        let mut string_losses = Vec::new();
        let mut endpoint_drag = Vec::new();
        let mut mass = identity(r);
        let mut bare_k = vec![0.0; r*r];
        for (i, b) in board.iter().enumerate() { bare_k[i*r+i] = (TAU*b.frequency_hz).powi(2); }
        let mut endpoint_k = bare_k;
        let mut loaded_add = vec![0.0; r*r];
        // Authored viscoelastic bending spectrum, not measured music wire.
        // Its loss is weighted by the actual bending/tension energy fraction.
        let bending = GeneralizedMaxwell::new(200e9, vec![(8e9, 0.0004), (2e9, 0.02)])
            .map_err(|e| e.to_string())?;
        let mut omitted_duplex_modes = 0;
        let mut capped_string_segments = 0;
        for (ci, c) in courses.iter().enumerate() {
            c.validate()?;
            for polarization in 0..planes {
            let group_start = strings.len();
            for member in 0..c.unison {
                let cents = (member as f64 - 0.5*(c.unison-1) as f64)*c.detune_cents;
                let tension = c.tension_at_cents(cents)?;
                for duplex in [false, true] {
                    if duplex && c.duplex_length_m == 0.0 { continue; }
                    let card = if duplex { Course { length_m: c.duplex_length_m, ..*c } } else { *c };
                    let retention = card.partial_retention(tension, band_hz, max_modes)?;
                    capped_string_segments += usize::from(retention.capped);
                    let start = modes.len();
                    let si = strings.len();
                    let mut beta2 = 0.0;
                    let mut modal_endpoint_k = 0.0;
                    let mut hammer_lift = c.strike_fraction;
                    let mut damper_lift = 0.35; // authored station, replace for measured dampers
                    let bridge: Vec<f64> = if polarization==0 {
                        board.iter().map(|b| b.bridge[usize::from(c.midi-21)]).collect()
                    } else {secondary.expect("admitted transverse rows")[ci].clone()};
                    let (mass_rate,strain_time)=if !damping {(0.0,0.0)}
                    else if rt0425_string_damping {
                        let (r_u,eta_u)=steinway_scale::string_damping_rt0425(c.midi)?;
                        (2.0*r_u,2.0*eta_u)
                    } else {(0.60,0.0)};
                    for n in 1..=retention.count {
                        let f = card.partial_hz(n, tension);
                        let omega = TAU*f;
                        let sign = if n % 2 == 0 { -1.0 } else { 1.0 };
                        let beta = det::sqrt(2.0*card.linear_density_kg_m*card.length_m)*sign/(n as f64*PI);
                        let hammer_shape = card.strike_shape(n);
                        let damper_shape = det::sin(n as f64*PI*0.35)/det::sqrt(card.modal_mass_kg());
                        let k = n as f64*PI/card.length_m;
                        let bend_fraction = card.flexural_rigidity_nm2*k*k/(tension + card.flexural_rigidity_nm2*k*k);
                        let strain_rate = if !damping {0.0} else if rt0425_string_damping {
                            strain_time*omega*omega
                        } else {omega*bending.loss_factor(omega)*bend_fraction};
                        string_losses.push(string_damping::ModeLoss {mass_rate,strain_rate});
                        oscillator.push(ModalAcousticMode { angular_frequency_rad_s: omega,
                            damping_ratio: 0.0, pressure_per_modal_velocity: C64::new(0.0,0.0) });
                        modes.push(StringMode { omega, beta, a: omega*omega*beta,
                            hammer_shape, damper_shape, string: si });
                        beta2 += beta*beta;
                        modal_endpoint_k += omega*omega*beta*beta;
                        hammer_lift -= hammer_shape*beta;
                        damper_lift -= damper_shape*beta;
                    }
                    if modes.len() == start {
                        if !duplex { return Err(format!("key {} fundamental exceeds the output band",c.midi)); }
                        omitted_duplex_modes += 1;
                    }
                    let residual_mass = card.linear_density_kg_m*card.length_m/3.0 - beta2;
                    if !residual_mass.is_finite() || residual_mass <= 0.0 { return Err("invalid residual string endpoint mass".into()); }
                    endpoint_drag.push(mass_rate*residual_mass+strain_time*tension/card.length_m);
                    for i in 0..r { for j in 0..r {
                        mass[i*r+j] += residual_mass*bridge[i]*bridge[j];
                        endpoint_k[i*r+j] += tension/card.length_m*bridge[i]*bridge[j];
                        loaded_add[i*r+j] += modal_endpoint_k*bridge[i]*bridge[j];
                    } }
                    let contact = if duplex || polarization!=0 { None } else {
                        let index = contact_strings.len(); contact_strings.push(si); Some(index)
                    };
                    strings.push(StringPort { course: ci, member, polarization, duplex, modes: start..modes.len(), bridge,
                        hammer_lift, damper_lift, contact });
                }
            }
            groups.push(group_start..strings.len());
            }
        }
        let loaded: Vec<f64> = endpoint_k.iter().zip(&loaded_add).map(|(x,y)| x+y).collect();
        let eig = fs_modal::eigh_gen_dense(&loaded, &mass, r).map_err(|e| e.to_string())?;
        for e in &eig {
            if !e.lambda.is_finite() || e.lambda <= 0.0 || e.residual > 1e-7*e.lambda
                || det::sqrt(e.lambda)/TAU > band_hz {
                return Err("loaded soundboard mode unresolved or above output band".into());
            }
        }
        let mut board_basis = vec![0.0; r*r];
        for i in 0..r { for j in 0..r { board_basis[i*r+j] = eig[j].phi[i]; } }
        for s in &mut strings {
            s.bridge = eig.iter().map(|e| s.bridge.iter().zip(&e.phi).map(|(g,p)| g*p).sum()).collect();
        }
        // Phi^T K Phi as two products, not one O(r^4) scalar expansion.
        let mut k_phi = vec![0.0; r*r];
        for a in 0..r { for j in 0..r {
            k_phi[a*r+j] = (0..r).map(|b| endpoint_k[a*r+b]*board_basis[b*r+j]).sum();
        } }
        let mut physical_board_k = vec![0.0; r*r];
        for i in 0..r { for j in i..r {
            let value = (0..r).map(|a| board_basis[a*r+i]*k_phi[a*r+j]).sum();
            physical_board_k[i*r+j] = value;
            physical_board_k[j*r+i] = value;
        } }
        let board_volume: Vec<f64> = eig.iter().map(|e| board.iter().zip(&e.phi).map(|(b,p)| b.volume*p).sum()).collect();
        let board_half_damping=if damping && board.iter().any(|b|b.damping_ratio>0.0) {
            Some(board_damping_flow(board,&board_basis,rate)?)
        } else {None};
        for e in &eig {
            // Board damping belongs to the physical matrix flow, not to an
            // invented mean loss on the loaded diagonal oscillator basis.
            oscillator.push(ModalAcousticMode { angular_frequency_rad_s: det::sqrt(e.lambda),
                damping_ratio: 0.0, pressure_per_modal_velocity: C64::new(0.0,0.0) });
        }
        let transition = transitions(rate, &oscillator)?;
        let n = modes.len();
        let string_loss=string_damping::Prepared::new(&strings,&modes,&groups,
            string_losses,endpoint_drag,r,0.5/f64::from(rate))?;
        let nc = contact_strings.len();
        let mut schur = vec![0.0; r*r];
        for i in 0..r { schur[i*r+i] = 1.0/transition[n+i].bq; }
        for s in &strings {
            let weight: f64 = s.modes.clone().map(|k| 0.25*modes[k].a*modes[k].a*transition[k].bq).sum();
            for i in 0..r { for j in 0..r { schur[i*r+j] -= weight*s.bridge[i]*s.bridge[j]; } }
        }
        let schur_inverse = inverse_spd(&schur,r)?;
        let mut contact_board = vec![0.0; nc*r];
        let mut contact_compliance = vec![0.0; nc*nc];
        for (c,&si) in contact_strings.iter().enumerate() {
            let s = &strings[si];
            let mut lift = s.hammer_lift;
            for k in s.modes.clone() {
                lift += 0.5*modes[k].a*transition[k].bq*modes[k].hammer_shape;
                contact_compliance[c*nc+c] += transition[k].bq*modes[k].hammer_shape.powi(2);
            }
            for j in 0..r { contact_board[c*r+j] = lift*s.bridge[j]; }
        }
        // W S^-1 W^T in O(nc*r^2 + nc^2*r), not O(nc^2*r^2).
        let mut response = vec![0.0; nc*r];
        for c in 0..nc { for a in 0..r {
            response[c*r+a] = (0..r).map(|b| schur_inverse[a*r+b]*contact_board[c*r+b]).sum();
        } }
        for i in 0..nc { for j in i..nc {
            let value: f64 = (0..r).map(|a| contact_board[i*r+a]*response[j*r+a]).sum();
            contact_compliance[i*nc+j] += value;
            if i != j { contact_compliance[j*nc+i] = contact_compliance[i*nc+j]; }
        } }
        Ok(Self { strings,groups,modes,contact_strings,board_count:r,q:vec![0.0;n+r],v:vec![0.0;n+r],
            next_q:vec![0.0;n+r],next_v:vec![0.0;n+r],contact_compliance,
            free_contact:vec![0.0;nc],rate,footprints:None,stretching:None,external_force:None,omitted_duplex_modes,capped_string_segments,transition,physical_board_k,
            board_half_damping,board_damping_enabled:damping,bare_board_damping_configured:false,
            string_loss,damped_v:vec![0.0;n+r],board_pre_loss:0.0,
            diagonal_omega2:oscillator.iter().map(|m|m.angular_frequency_rad_s.powi(2)).collect(),last_modal_loss_j:0.0,
            board_volume,board_basis,schur_inverse,contact_response:response,
            #[cfg(test)]
            contact_board,
            free_q:vec![0.0;n+r],free_v:vec![0.0;n+r],
            r_string:vec![0.0;n],board_rhs:vec![0.0;r],board_end:vec![0.0;r] })
    }

    /// Replace bare-board modal loss with a complete supplied viscous matrix
    /// [1/s], row-major in the original mass-normalized bare-board coordinates.
    /// A Ritz basis generally needs this dense form even when source damping
    /// was diagonal. It is not added to BoardMode::damping_ratio.
    ///
    /// Cold and once only. Input must be finite, exactly symmetric and resolved
    /// PSD; the same scale-aware nullspace tolerance as the existing damping
    /// flow applies. Return the effective continuous Phi^T C Phi used by that
    /// flow, so harmonic loss uses the same physical operator. Global damping
    /// false still validates the supplied matrix but returns an all-zero image.
    pub fn configure_bare_board_damping(&mut self,c:&[f64])->Result<Vec<f64>,String> {
        let r=self.board_count;
        if self.bare_board_damping_configured || self.q.iter().chain(&self.v).any(|x|*x!=0.0) {
            return Err("bare-board damping requires a fresh unconfigured bank".into());
        }
        if c.len()!=r*r || c.iter().any(|x|!x.is_finite())
            || (0..r).any(|i|(i+1..r).any(|j|c[i*r+j]!=c[j*r+i])) {
            return Err("bare-board damping must be a complete finite symmetric matrix".into());
        }
        // Validate the physical input too: an ill-scaled congruence must not
        // conceal a negative source direction, including in lossless mode.
        board_damping_spectrum(c,r)?;
        let mut c_phi=vec![0.0;r*r];
        for a in 0..r {for j in 0..r {
            c_phi[a*r+j]=(0..r).map(|b|c[a*r+b]*self.board_basis[b*r+j]).sum();
        }}
        let mut loaded=vec![0.0;r*r];
        for i in 0..r {for j in i..r {
            let value=(0..r).map(|a|self.board_basis[a*r+i]*c_phi[a*r+j]).sum();
            loaded[i*r+j]=value;loaded[j*r+i]=value;
        }}
        let flow=board_damping_matrix_flow(&loaded,r,self.rate)?;
        // Publication follows complete admission. The geometry, conservative
        // map, contact compliance and all authoritative motion remain intact.
        self.board_half_damping=if self.board_damping_enabled && loaded.iter().any(|x|*x!=0.0) {
            Some(flow)
        } else {None};
        if !self.board_damping_enabled {loaded.fill(0.0);}
        self.bare_board_damping_configured=true;
        Ok(loaded)
    }

    pub fn has_secondary_polarization(&self)->bool {
        self.strings.iter().any(|s|s.polarization==1)
    }

    /// Cold projection of a physical bare-board shape into the SAME loaded
    /// coordinates used for bridge work. This preserves modal normalization for
    /// spatial microphones; a unit-mass or frequency-only observer is not used.
    pub fn project_board_shape(&self, shape: &[f64]) -> Result<Vec<f64>, String> {
        let r=self.board_count;
        if shape.len()!=r || shape.iter().any(|x| !x.is_finite()) {
            return Err("surface projection must cover every finite bare-board mode".into());
        }
        let out:Vec<f64>=(0..r).map(|j|(0..r).map(|i|shape[i]*self.board_basis[i*r+j]).sum()).collect();
        if out.iter().any(|x| !x.is_finite()) { return Err("surface projection overflow".into()); }
        Ok(out)
    }

    /// The exact row-major map used by `project_board_shape`, followed by the
    /// loaded board's diagonal reference omega-squared values. Rows are bare
    /// coordinates; columns are the pressure trace's loaded coordinates.
    /// These are split-step reference oscillators, not full instrument poles.
    pub fn board_pressure_basis(&self) -> Result<(&[f64], &[f64]), String> {
        let r = self.board_count;
        let omega2 = self.diagonal_omega2.get(self.modes.len()..)
            .ok_or("loaded board reference frequencies are absent")?;
        if !(1..=MAX_BOARD_MODES).contains(&r) || self.board_basis.len() != r*r
            || omega2.len() != r || self.board_basis.iter().any(|x| !x.is_finite())
            || omega2.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err("loaded pressure basis dimensions or values are invalid".into());
        }
        Ok((&self.board_basis, omega2))
    }

    pub fn contact_position(&self, contact: usize, q: &[f64]) -> f64 {
        let s = &self.strings[self.contact_strings[contact]];
        let n = self.modes.len();
        if let Some(p) = &self.footprints {
            let p = &p.points[contact];
            let primary = s.modes.clone().zip(&p.shape).map(|(k,g)| g*q[k]).sum::<f64>()
                + p.lift*s.bridge.iter().zip(&q[n..]).map(|(g,x)|g*x).sum::<f64>();
            let mut displacement = p.normal[0]*primary;
            if let Some(side) = p.secondary {
                let s = &self.strings[side];
                let lateral = s.modes.clone().zip(&p.shape).map(|(k,g)| g*q[k]).sum::<f64>()
                    + p.lift*s.bridge.iter().zip(&q[n..]).map(|(g,x)|g*x).sum::<f64>();
                displacement += p.normal[1]*lateral;
            }
            return displacement;
        }
        s.modes.clone().map(|k| self.modes[k].hammer_shape*q[k]).sum::<f64>()
            + s.hammer_lift*s.bridge.iter().zip(&q[n..]).map(|(g,x)| g*x).sum::<f64>()
    }

    /// Cold admission of additional force ports. The force is numerical trial
    /// scratch; its owner retains physical contact and material history.
    pub(super) fn prepare_external_force(&mut self) {
        self.external_force = Some(vec![0.0; self.q.len()]);
    }

    pub(super) fn set_external_force(&mut self, force: &[f64]) -> Result<(), &'static str> {
        let prepared = self.external_force.as_mut().ok_or("external force ports were not prepared")?;
        if prepared.len() != force.len() || force.iter().any(|x| !x.is_finite()) {
            return Err("external force must cover every finite bank coordinate");
        }
        prepared.copy_from_slice(force);
        Ok(())
    }

    /// Displacement response to a held generalized load in this exact ZOH/
    /// Schur image. Cold contact preparation uses g_i^T H g_j, including the
    /// moving bridge and sympathetic strings. No state or clock advances.
    ///
    /// Eliminating the string responses gives
    /// S db = f_b + sum(g * a * bq * f_s / 2), then
    /// ds = bq * (f_s + a * g^T db / 2). This is the same symmetric H used by
    /// predict/finish, not an independent contact compliance approximation.
    pub(super) fn force_response(&self, force: &[f64], out: &mut [f64]) -> Result<(), String> {
        let n = self.modes.len();
        let r = self.board_count;
        if force.len() != n+r || out.len() != n+r || force.iter().any(|x| !x.is_finite()) {
            return Err("force response needs the complete finite bank coordinates".into());
        }
        let mut rhs = force[n..].to_vec();
        for group in &self.groups {
            let g = &self.strings[group.start].bridge;
            let mut reaction = 0.0;
            for si in group.clone() { for k in self.strings[si].modes.clone() {
                reaction += 0.5*self.modes[k].a*self.transition[k].bq*force[k];
            }}
            for (value, shape) in rhs.iter_mut().zip(g) { *value += shape*reaction; }
        }
        for j in 0..r {
            out[n+j] = (0..r).map(|k| self.schur_inverse[j*r+k]*rhs[k]).sum();
        }
        for group in &self.groups {
            let g = &self.strings[group.start].bridge;
            let bridge = g.iter().zip(&out[n..]).map(|(g,q)| g*q).sum::<f64>();
            for si in group.clone() { for k in self.strings[si].modes.clone() {
                out[k] = self.transition[k].bq*(force[k]+0.5*self.modes[k].a*bridge);
            }}
        }
        if out.iter().any(|x| !x.is_finite()) { return Err("force response overflow".into()); }
        Ok(())
    }

    /// Prepare the unforced coupled end positions. No allocations.
    pub fn predict(&mut self) {
        if let Some(stretching)=&mut self.stretching {stretching.predict();}
        let n = self.modes.len(); let r = self.board_count;
        self.damped_v.copy_from_slice(&self.v);
        self.board_pre_loss=self.string_loss.apply(&self.strings,&self.modes,
            &self.groups,&mut self.damped_v);
        if let Some(flow)=&self.board_half_damping {
            self.board_pre_loss+=damp_board(flow,&self.damped_v[n..],&mut self.board_rhs);
            self.damped_v[n..].copy_from_slice(&self.board_rhs);
        }
        for k in 0..self.q.len() {
            let t=self.transition[k];
            let v=self.damped_v[k];
            self.free_q[k]=t.qq*self.q[k]+t.qv*v;
            self.free_v[k]=t.vq*self.q[k]+t.vv*v;
            if let Some(stretching)=&self.stretching {
                // Held trial force; updated together with contact before ANY
                // state publication. The base exact-ZOH map is unchanged.
                self.free_q[k]+=t.bq*stretching.force[k];
                self.free_v[k]+=t.bv*stretching.force[k];
            }
            if let Some(force)=&self.external_force {
                self.free_q[k]+=t.bq*force[k];
                self.free_v[k]+=t.bv*force[k];
            }
        }
        for j in 0..r { self.board_rhs[j]=self.free_q[n+j]/self.transition[n+j].bq; }
        for group in &self.groups {
            let g=&self.strings[group.start].bridge;
            let b0=g.iter().zip(&self.q[n..]).map(|(g,q)| g*q).sum::<f64>();
            let mut reaction=0.0;
            for si in group.clone() { for k in self.strings[si].modes.clone() {
                let m=self.modes[k];
                let rs=self.free_q[k]+0.5*self.transition[k].bq*m.a*b0;
                self.r_string[k]=rs;
                reaction+=0.5*m.a*(self.q[k]+rs);
            } }
            for (rhs,g) in self.board_rhs.iter_mut().zip(g) { *rhs+=g*reaction; }
        }
        for j in 0..r {
            self.board_end[j]=(0..r).map(|k| self.schur_inverse[j*r+k]*self.board_rhs[k]).sum();
            self.next_q[n+j]=self.board_end[j];
        }
        for group in &self.groups {
            let g=&self.strings[group.start].bridge;
            let b1=g.iter().zip(&self.board_end).map(|(g,q)| g*q).sum::<f64>();
            for si in group.clone() {
                let s=&self.strings[si];
                let mut position=s.hammer_lift*b1;
                for k in s.modes.clone() {
                    let m=self.modes[k];
                    self.next_q[k]=self.r_string[k]+0.5*self.transition[k].bq*m.a*b1;
                    position+=m.hammer_shape*self.next_q[k];
                }
                if let Some(c)=s.contact { self.free_contact[c]=position; }
            }
        }
        if self.footprints.is_some() {
            for c in 0..self.free_contact.len() {
                self.free_contact[c] = self.contact_position(c, &self.next_q);
            }
        }
    }

    /// Finish the SAME prediction with held contact forces [N]. Candidate
    /// arrays are not committed until the nonlinear island accepts the step.
    pub fn finish(&mut self, forces: &[f64]) {
        let n=self.modes.len();let r=self.board_count;
        self.last_modal_loss_j=self.board_pre_loss;
        if let Some(stretching)=&self.stretching {self.last_modal_loss_j+=stretching.modal_loss();}
        assert_eq!(forces.len(),self.contact_strings.len(),"one force per contact");
        if let Some(p) = &mut self.footprints { p.gather(&self.strings, forces); }
        // S^-1 W^T was already computed for the compliance matrix. Keep it
        // rather than repeat the small dense solve for every accepted substep.
        // Silent contacts skip only multiplication, not sympathetic dynamics.
        self.next_q[n..].copy_from_slice(&self.board_end);
        for (c,&force) in forces.iter().enumerate() {
            if force==0.0 { continue; }
            for j in 0..r { self.next_q[n+j]+=self.contact_response[c*r+j]*force; }
        }
        self.board_rhs.fill(0.0);
        for group in &self.groups {
            let g=&self.strings[group.start].bridge;
            let bbar=g.iter().enumerate().map(|(j,g)|g*0.5*(self.q[n+j]+self.next_q[n+j])).sum::<f64>();
            let mut reaction=0.0;
            for si in group.clone() {
                let s=&self.strings[si];
                let contact=s.contact.map_or(0.0,|c|forces[c]);
                reaction+=self.footprints.as_ref().map_or(s.hammer_lift*contact, |p|p.bridge_force[si]);
                for k in s.modes.clone() {
                    let m=self.modes[k];
                    let f=m.a*bbar+self.footprints.as_ref().map_or(contact*m.hammer_shape, |p|p.modal_force[k]);
                    self.next_q[k]=self.free_q[k]+self.transition[k].bq*f;
                    self.next_v[k]=self.free_v[k]+self.transition[k].bv*f;
                    reaction+=m.a*0.5*(self.q[k]+self.next_q[k]);
                    let mut total_f=self.stretching.as_ref().map_or(f,|s|f+s.force[k]);
                    if let Some(force)=&self.external_force {total_f+=force[k];}
                    self.last_modal_loss_j+=total_f*(self.next_q[k]-self.q[k])-0.5*(
                        self.next_v[k].powi(2)-self.damped_v[k].powi(2)
                        +self.diagonal_omega2[k]*(self.next_q[k].powi(2)-self.q[k].powi(2)));
                }
            }
            for (rhs,g) in self.board_rhs.iter_mut().zip(g) { *rhs+=g*reaction; }
        }
        for j in 0..r {
            let k=n+j;let f=self.board_rhs[j];
            self.next_v[k]=self.free_v[k]+self.transition[k].bv*f;
            let mut total_f=self.stretching.as_ref().map_or(f,|s|f+s.force[k]);
            if let Some(force)=&self.external_force {total_f+=force[k];}
            self.last_modal_loss_j+=total_f*(self.next_q[k]-self.q[k])-0.5*(
                self.next_v[k].powi(2)-self.damped_v[k].powi(2)
                +self.diagonal_omega2[k]*(self.next_q[k].powi(2)-self.q[k].powi(2)));
        }
        if let Some(flow)=&self.board_half_damping {
            self.last_modal_loss_j+=damp_board(flow,&self.next_v[n..],&mut self.board_rhs);
            self.next_v[n..].copy_from_slice(&self.board_rhs);
        }
        self.last_modal_loss_j+=self.string_loss.apply(&self.strings,&self.modes,
            &self.groups,&mut self.next_v);
    }

    pub fn commit(&mut self) {
        self.q.copy_from_slice(&self.next_q);self.v.copy_from_slice(&self.next_v);
        if let Some(stretching)=&mut self.stretching {stretching.commit();}
    }

    pub fn energy(&self) -> f64 { self.energy_at(&self.q,&self.v) }
    /// Evaluate transverse/board coordinates with the authoritative axial
    /// state held fixed. Use next_energy() for the complete uncommitted trial.
    pub fn energy_at(&self,q:&[f64],v:&[f64])->f64 {
        let mut energy=self.linear_energy_at(q,v);
        if let Some(s)=&self.stretching {energy+=s.energy(q,&self.strings,&self.modes);}
        energy
    }
    pub fn next_energy(&self)->f64 {
        let mut energy=self.linear_energy_at(&self.next_q,&self.next_v);
        if let Some(s)=&self.stretching {energy+=s.candidate_energy(&self.next_q,&self.strings,&self.modes);}
        energy
    }
    fn linear_energy_at(&self,q:&[f64],v:&[f64])->f64 {
        let n=self.modes.len();let r=self.board_count;
        let mut energy=0.5*v.iter().map(|x|x*x).sum::<f64>();
        for group in &self.groups {
            let b=self.strings[group.start].bridge.iter().zip(&q[n..]).map(|(g,q)|g*q).sum::<f64>();
            for si in group.clone() { for k in self.strings[si].modes.clone() {
                let m=self.modes[k];
                energy+=0.5*m.omega.powi(2)*(q[k]-m.beta*b).powi(2);
            } }
        }
        for i in 0..r { for j in 0..r {energy+=0.5*q[n+i]*self.physical_board_k[i*r+j]*q[n+j];} }
        energy
    }

    /// Everything a contact-free mechanics step of this bank reads, as plain data, for an external
    /// stepper (the GPU renderer in apps/instrument-lab/render-gpu-piano). Exporting changes nothing; the
    /// bank's own predict / finish / commit stay the reference. None when an option the external stepper
    /// does not model (string stretching, hammer footprints, held external forces) is active.
    #[allow(dead_code)] // used by apps/instrument-lab/render-gpu-piano
    pub fn gpu_image(&self)->Option<BankImage> {
        if self.stretching.is_some() || self.footprints.is_some() || self.external_force.is_some() { return None; }
        Some(BankImage {
            n: self.modes.len(), r: self.board_count,
            transition: self.transition.iter().map(|t| [t.qq, t.qv, t.vq, t.vv, t.bq, t.bv]).collect(),
            modes: self.modes.clone(), omega2: self.diagonal_omega2.clone(),
            strings: self.strings.clone(), groups: self.groups.clone(),
            schur_inverse: self.schur_inverse.clone(),
            board_half_damping: self.board_half_damping.clone(),
            string_loss: self.string_loss.image(),
            contact_strings: self.contact_strings.clone(),
            q: self.q.clone(), v: self.v.clone(),
        })
    }

    /// A free-oscillation estimate bounding |contact_position(contact)|: each string mode's and board
    /// coordinate's amplitude sqrt(q^2+(v/omega)^2), weighted by the contact's shape. A planner's
    /// estimate only; the coupled, damped motion is not bounded by it exactly.
    #[allow(dead_code)] // used by apps/instrument-lab/render-gpu-piano
    pub fn contact_amplitude_bound(&self, contact: usize) -> f64 {
        let s=&self.strings[self.contact_strings[contact]];let n=self.modes.len();
        let amplitude=|k:usize| det::sqrt(self.q[k]*self.q[k]+self.v[k]*self.v[k]/self.diagonal_omega2[k]);
        s.modes.clone().map(|k| self.modes[k].hammer_shape.abs()*amplitude(k)).sum::<f64>()
            + s.hammer_lift.abs()*s.bridge.iter().enumerate().map(|(j,g)| g.abs()*amplitude(n+j)).sum::<f64>()
    }

    pub fn volume_velocity(&self)->f64 {
        self.board_volume.iter().zip(&self.v[self.modes.len()..]).map(|(a,v)|a*v).sum()
    }

    /// Exact dissipative rank-one velocity map at an actual displacement port.
    /// This is a viscous damper approximation, not a measured wool-pad model.
    /// `drag` is N s/m; damper/rank-one flow removes kinetic energy only.
    pub fn damp_string(&mut self,si:usize,drag:f64,dt:f64)->f64 {
        let s=&self.strings[si];let modes=&self.modes;
        dampers::rank_one(&mut self.v,s.modes.clone(),modes.len(),
            |k|modes[k].damper_shape,&s.bridge,s.damper_lift,drag,dt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn published_string_loss_is_used_by_the_prepared_bank_without_changing_frequencies() {
        let c = super::super::steinway_scale::courses().unwrap()[48]; // A4
        let board = super::super::board::demonstration();
        let old = Bank::new(&[c], &board, 192_000, 21_600.0, 4, true).unwrap();
        let source = Bank::new_with_string_damping(&[c], &board, 192_000,
            21_600.0, 4, true, None, true).unwrap();
        assert_eq!(old.modes[0].omega, source.modes[0].omega);
        let omega = source.modes[0].omega;
        let (r_u, eta_u) = steinway_scale::string_damping_rt0425(69).unwrap();
        let loss=source.string_loss.modes[0];
        assert_eq!(loss.mass_rate,2.0*r_u);
        assert_eq!(loss.strain_rate,2.0*eta_u*omega*omega);
        assert!(((loss.mass_rate+loss.strain_rate)/(2.0*omega)
            -(r_u/omega+eta_u*omega)).abs()<1e-14);
        assert_eq!(source.transition[0].vv,old.transition[0].vv);
        assert_ne!(loss.strain_rate,old.string_loss.modes[0].strain_rate);
    }
    fn bank(damped:bool)->Bank{
        let scale=super::super::geometry::demonstration_scale().unwrap();
        Bank::new(&[scale[48]],&super::super::board::demonstration(),192_000,21_600.0,12,damped).unwrap()
    }
    #[test]
    fn reciprocal_bridge_conserves_energy_without_gyroscopic_coupling(){
        let mut b=bank(false);
        for k in 0..b.q.len(){b.q[k]=1e-6*det::sin(k as f64);b.v[k]=0.01*det::cos(k as f64);}
        let initial=b.energy();let forces=vec![0.0;b.contact_strings.len()];
        for _ in 0..1000{b.predict();b.finish(&forces);b.commit();}
        assert!((b.energy()/initial-1.0).abs()<1e-8);
    }
    #[test]
    fn contact_compliance_is_reciprocal_and_force_work_is_physical(){
        let mut b=bank(true);let nc=b.contact_strings.len();
        for i in 0..nc{assert!(b.contact_compliance[i*nc+i]>0.0);for j in 0..nc{
            assert!((b.contact_compliance[i*nc+j]-b.contact_compliance[j*nc+i]).abs()<1e-15);
        }}
        let forces:Vec<f64>=(1..=nc).map(|x|x as f64).collect();
        let before=b.energy();b.predict();b.finish(&forces);
        let work=(0..nc).map(|i|forces[i]*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum::<f64>();
        let defect=b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work;
        assert!(defect.abs()<1e-12,"{defect:e}");
    }
    #[test]
    fn repeated_contact_trials_restart_from_the_same_physically_damped_state() {
        let mut b=bank(true);
        for k in 0..b.q.len() {b.q[k]=1e-8*det::sin(k as f64);b.v[k]=1e-3*det::cos(k as f64);}
        let q=b.q.clone();let v=b.v.clone();let before=b.energy();
        let forces=vec![0.05;b.contact_strings.len()];
        b.predict();b.finish(&forces);
        let accepted_q=b.next_q.clone();let accepted_v=b.next_v.clone();let loss=b.last_modal_loss_j;
        for predict_again in [false,true] {
            if predict_again {b.predict();}
            b.finish(&forces);
            assert_eq!(b.q,q);assert_eq!(b.v,v);
            assert_eq!(b.next_q,accepted_q);assert_eq!(b.next_v,accepted_v);
            assert_eq!(b.last_modal_loss_j,loss);
        }
        let work: f64=forces.iter().enumerate().map(|(i,f)|
            f*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum();
        assert!(loss>0.0);
        assert!((b.next_energy()-before+loss-work).abs()<1e-12);
    }
    #[test]
    fn generalized_force_response_is_reciprocal_and_matches_contact_ports() {
        let b=bank(true);let d=b.q.len();let nc=b.contact_strings.len();
        let a:Vec<_>=(0..d).map(|k|det::sin(0.31*k as f64)).collect();
        let c:Vec<_>=(0..d).map(|k|det::cos(0.73*k as f64)).collect();
        let mut ha=vec![0.0;d];let mut hc=vec![0.0;d];
        b.force_response(&a,&mut ha).unwrap();b.force_response(&c,&mut hc).unwrap();
        let ahc:f64=a.iter().zip(&hc).map(|(a,h)|a*h).sum();
        let cha:f64=c.iter().zip(&ha).map(|(c,h)|c*h).sum();
        assert!((ahc-cha).abs()<1e-20+1e-12*ahc.abs());
        assert!(a.iter().zip(&ha).map(|(a,h)|a*h).sum::<f64>()>0.0);
        // Independently extract physical point-displacement rows and compare
        // against the pre-existing hammer contact compliance owner.
        let mut unit=vec![0.0;d];let mut g=vec![0.0;d];
        for j in 0..nc {
            for k in 0..d {
                unit[k]=1.0;g[k]=b.contact_position(j,&unit);unit[k]=0.0;
            }
            b.force_response(&g,&mut ha).unwrap();
            for i in 0..nc {
                assert!((b.contact_position(i,&ha)-b.contact_compliance[i*nc+j]).abs()<1e-15);
            }
        }
        assert!(b.force_response(&a[..d-1],&mut ha).is_err());
        assert!(b.force_response(&vec![f64::NAN;d],&mut ha).is_err());
        assert!(b.q.iter().chain(&b.v).all(|x|*x==0.0));
    }
    #[test]
    fn generalized_loads_share_the_hammer_step_and_close_complete_force_work() {
        for damped in [false,true] {
            let mut b=bank(damped);let d=b.q.len();let nc=b.contact_strings.len();
            for k in 0..d {b.q[k]=1e-8*det::sin(k as f64);b.v[k]=1e-4*det::cos(k as f64);}
            let hammer:Vec<_>=(0..nc).map(|i|0.1*(i+1) as f64).collect();
            let mut force=vec![0.0;d];let mut response=vec![0.0;d];
            assert!(b.set_external_force(&force).is_err());
            b.prepare_external_force();
            for step in 0..48 {
                force.fill(0.0);b.set_external_force(&force).unwrap();
                b.predict();b.finish(&hammer);let free=b.next_q.clone();
                for (k,f) in force.iter_mut().enumerate() {*f=0.02*det::sin((k+3*step) as f64*0.31);}
                b.force_response(&force,&mut response).unwrap();
                b.set_external_force(&force).unwrap();
                assert!(b.set_external_force(&[f64::NAN]).is_err());
                assert_eq!(b.external_force.as_deref(),Some(force.as_slice()));
                let before=b.energy();b.predict();b.finish(&hammer);
                for k in 0..d {
                    assert!((b.next_q[k]-free[k]-response[k]).abs()<1e-16);
                }
                let work:f64=force.iter().enumerate().map(|(k,f)|f*(b.next_q[k]-b.q[k])).sum::<f64>()
                    +hammer.iter().enumerate().map(|(i,f)|f*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum::<f64>();
                let defect=b.next_energy()-before+b.last_modal_loss_j-work;
                assert!(defect.abs()<1e-12,"held-load work defect {defect:e}");
                assert!(b.last_modal_loss_j>=-1e-13);b.commit();
            }
        }
    }
    #[test]
    fn damper_port_removes_exactly_its_reported_energy(){
        let mut b=bank(false);b.v.fill(0.03);let before=b.energy();
        let loss=b.damp_string(0,0.4,1.0/192_000.0);
        assert!(loss>0.0);assert!((b.energy()+loss-before).abs()<1e-12);
    }
    #[test]
    fn surface_projection_uses_the_same_loaded_basis_as_volume_velocity(){
        let mut b=bank(false);
        let shape:Vec<f64>=super::super::board::demonstration().iter().map(|m|m.volume).collect();
        let projected=b.project_board_shape(&shape).unwrap();
        assert_eq!(projected,b.board_volume);
        b.v.fill(0.03);
        let observed=projected.iter().zip(&b.v[b.modes.len()..]).map(|(g,v)|g*v).sum::<f64>();
        assert_eq!(observed,b.volume_velocity());
        assert!(b.project_board_shape(&[]).is_err());
    }

    #[test]
    fn exported_pressure_basis_reproduces_projection_for_register_and_contact_arms() {
        let scale = super::super::geometry::demonstration_scale().unwrap();
        let board = super::super::board::demonstration();
        for key in [60_u8, 69, 84] {
            let course = *scale.iter().find(|c| c.midi == key).unwrap();
            for secondary in [None, Some(vec![vec![0.4; board.len()]])] {
                let b = Bank::new_with_transverse_bridge(&[course], &board,
                    192_000, 21_600.0, 12, true, secondary.as_deref()).unwrap();
                let (map, omega2) = b.board_pressure_basis().unwrap();
                for i in 0..board.len() {
                    let mut shape = vec![0.0; board.len()]; shape[i] = 1.0;
                    assert_eq!(b.project_board_shape(&shape).unwrap(),
                        map[i*board.len()..(i+1)*board.len()]);
                }
                assert_eq!(omega2, &b.diagonal_omega2[b.modes.len()..]);
                assert_ne!(omega2[0], (TAU*board[0].frequency_hz).powi(2));
            }
        }
    }

    #[test]
    fn pressure_basis_refuses_missing_mismatched_and_nonfinite_coordinates() {
        let mut b = bank(false);
        let original = b.board_basis.clone();
        b.board_basis.pop(); assert!(b.board_pressure_basis().is_err());
        b.board_basis = original; b.board_basis[0] = f64::NAN;
        assert!(b.board_pressure_basis().is_err());
        b.board_basis[0] = 0.0; b.diagonal_omega2.pop();
        assert!(b.board_pressure_basis().is_err());
        b.diagonal_omega2.push(-1.0); assert!(b.board_pressure_basis().is_err());
    }

    /// The former scalar expansion, independent of the course grouping.
    /// Kept only as a small direct regression oracle, not another runtime image.
    fn unfactored(b:&Bank, forces:&[f64])->(Vec<f64>,Vec<f64>,Vec<f64>,f64) {
        let n=b.modes.len();let r=b.board_count;
        let mut v0=b.v.clone();let mut loss=0.0;
        let mut string_loss=b.string_loss.clone();
        loss+=string_loss.apply(&b.strings,&b.modes,&b.groups,&mut v0);
        if let Some(flow)=&b.board_half_damping {
            let before=v0[n..].to_vec();loss+=damp_board(flow,&before,&mut v0[n..]);
        }
        let fq:Vec<f64>=(0..n+r).map(|k|b.transition[k].qq*b.q[k]+b.transition[k].qv*v0[k]).collect();
        let fv:Vec<f64>=(0..n+r).map(|k|b.transition[k].vq*b.q[k]+b.transition[k].vv*v0[k]).collect();
        let mut rhs:Vec<f64>=(0..r).map(|j|fq[n+j]/b.transition[n+j].bq).collect();
        let mut rs=vec![0.0;n];
        for (k,m) in b.modes.iter().enumerate() {
            let g=&b.strings[m.string].bridge;
            let b0=g.iter().zip(&b.q[n..]).map(|(g,q)| g*q).sum::<f64>();
            rs[k]=fq[k]+0.5*b.transition[k].bq*m.a*b0;
            for j in 0..r {rhs[j]+=0.5*m.a*g[j]*(b.q[k]+rs[k]);}
        }
        let end:Vec<f64>=(0..r).map(|j|(0..r).map(|a|b.schur_inverse[j*r+a]*rhs[a]).sum()).collect();
        let mut q=vec![0.0;n+r];let mut v=vec![0.0;n+r];q[n..].copy_from_slice(&end);
        for (k,m) in b.modes.iter().enumerate() {
            let b1=b.strings[m.string].bridge.iter().zip(&end).map(|(g,q)| g*q).sum::<f64>();
            q[k]=rs[k]+0.5*b.transition[k].bq*m.a*b1;
        }
        let free=(0..forces.len()).map(|c|b.contact_position(c,&q)).collect();
        for j in 0..r {for a in 0..r {
            let f=(0..forces.len()).map(|c|b.contact_board[c*r+a]*forces[c]).sum::<f64>();
            q[n+j]+=b.schur_inverse[j*r+a]*f;
        }}
        for (k,m) in b.modes.iter().enumerate() {
            let s=&b.strings[m.string];
            let bbar=(0..r).map(|j|s.bridge[j]*0.5*(b.q[n+j]+q[n+j])).sum::<f64>();
            let f=m.a*bbar+s.contact.map_or(0.0,|c|forces[c]*m.hammer_shape);
            q[k]=fq[k]+b.transition[k].bq*f;v[k]=fv[k]+b.transition[k].bv*f;
            loss+=f*(q[k]-b.q[k])-0.5*(v[k]*v[k]-v0[k]*v0[k]+b.diagonal_omega2[k]*(q[k]*q[k]-b.q[k]*b.q[k]));
        }
        rhs.fill(0.0);
        for (k,m) in b.modes.iter().enumerate() {for j in 0..r {
            rhs[j]+=m.a*b.strings[m.string].bridge[j]*0.5*(b.q[k]+q[k]);
        }}
        for (c,&si) in b.contact_strings.iter().enumerate() {for j in 0..r {
            rhs[j]+=b.strings[si].hammer_lift*b.strings[si].bridge[j]*forces[c];
        }}
        for j in 0..r {
            let k=n+j;v[k]=fv[k]+b.transition[k].bv*rhs[j];
            loss+=rhs[j]*(q[k]-b.q[k])-0.5*(v[k]*v[k]-v0[k]*v0[k]+b.diagonal_omega2[k]*(q[k]*q[k]-b.q[k]*b.q[k]));
        }
        if let Some(flow)=&b.board_half_damping {
            let original=v[n..].to_vec();loss+=damp_board(flow,&original,&mut v[n..]);
        }
        loss+=string_loss.apply(&b.strings,&b.modes,&b.groups,&mut v);
        (free,q,v,loss)
    }
    fn close(a:&[f64],b:&[f64]) {
        assert_eq!(a.len(),b.len());
        for (x,y) in a.iter().zip(b) {assert!((x-y).abs()<2e-12*(1.0+x.abs().max(y.abs())),"{x:e} != {y:e}");}
    }
    #[test]
    fn grouped_bridge_matches_unfactored_dynamics_with_unisons_and_duplexes() {
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let courses=[scale[0],scale[39],scale[48]];
        for damped in [false,true] {
            let mut b=Bank::new(&courses,&super::super::board::demonstration(),192_000,21_600.0,24,damped).unwrap();
            assert_eq!(b.groups.len(),courses.len());
            for group in &b.groups {for si in group.clone() {
                assert_eq!(b.strings[si].bridge,b.strings[group.start].bridge);
            }}
            for k in 0..b.q.len(){b.q[k]=1e-9*det::sin(k as f64);b.v[k]=1e-4*det::cos(k as f64);}
            for sample in 0..64 {
                let force:Vec<f64>=(0..b.contact_strings.len()).map(|c|0.5*det::sin((c+sample) as f64)).collect();
                let (free,q,v,loss)=unfactored(&b,&force);
                b.predict();close(&free,&b.free_contact);b.finish(&force);
                close(&q,&b.next_q);close(&v,&b.next_v);
                assert!((loss-b.last_modal_loss_j).abs()<1e-12);
                b.commit();
            }
        }
    }
    #[test]
    fn factored_contact_operator_reproduces_each_unit_force_column() {
        let mut b=bank(false);let nc=b.contact_strings.len();let mut forces=vec![0.0;nc];
        for c in 0..nc {
            forces.fill(0.0);forces[c]=1.0;b.predict();b.finish(&forces);
            for i in 0..nc {
                let actual=b.contact_position(i,&b.next_q)-b.free_contact[i];
                assert!((actual-b.contact_compliance[i*nc+c]).abs()<1e-14);
            }
        }
    }

    #[test]
    fn physical_board_damping_preserves_heterogeneous_modal_losses_and_mass_loading() {
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let mut board=super::super::board::demonstration();
        for (m,zeta) in board.iter_mut().zip([0.0,0.001,0.07,0.2]) {m.damping_ratio=zeta;}
        let b=Bank::new(&[scale[48]],&board,192_000,21_600.0,24,true).unwrap();
        let r=board.len();let v:Vec<f64>=(0..r).map(|i|0.03*(i+1) as f64).collect();
        let expected_power: f64=board.iter().enumerate().map(|(i,m)| {
            let speed: f64=(0..r).map(|j|b.board_basis[i*r+j]*v[j]).sum();
            2.0*m.damping_ratio*TAU*m.frequency_hz*speed*speed
        }).sum();
        let mut out=vec![0.0;r];
        let loss=damp_board(b.board_half_damping.as_ref().unwrap(),&v,&mut out);
        assert!(loss>0.0);
        // Compare the exact finite-time loss to its physical infinitesimal power.
        let ratio=loss/(0.5/f64::from(b.rate)*expected_power);
        assert!((0.99..=1.000001).contains(&ratio),"{ratio}");
        assert!((loss+0.5*out.iter().map(|x|x*x).sum::<f64>()
            -0.5*v.iter().map(|x|x*x).sum::<f64>()).abs()<1e-15);
    }
    #[test]
    fn dense_bare_board_damping_keeps_rank_one_flow_and_complete_force_work() {
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let board=super::super::board::demonstration();let r=board.len();assert_eq!(r,4);
        let mut b=Bank::new(&[scale[48]],&board,192_000,21_600.0,12,true).unwrap();
        let stiffness=b.physical_board_k.clone();let compliance=b.contact_compliance.clone();
        let basis=b.board_basis.clone();
        let u=[4.0,-3.0,2.0,1.0];let mut c=vec![0.0;r*r];
        for i in 0..r {for j in 0..r {c[i*r+j]=2.0*u[i]*u[j];}}
        let loaded=b.configure_bare_board_damping(&c).unwrap();
        assert_eq!(b.physical_board_k,stiffness);assert_eq!(b.contact_compliance,compliance);
        assert_eq!(b.board_basis,basis);assert!(b.q.iter().chain(&b.v).all(|x|*x==0.0));
        let axis=b.project_board_shape(&u).unwrap();
        for i in 0..r {for j in 0..r {
            assert!((loaded[i*r+j]-2.0*axis[i]*axis[j]).abs()<1e-12);
        }}
        let v=[0.02,-0.013,0.007,0.025];let mut out=vec![0.0;r];
        let norm2=axis.iter().map(|x|x*x).sum::<f64>();
        let along=axis.iter().zip(&v).map(|(a,b)|a*b).sum::<f64>()/norm2;
        // Independent closed form for exp(-dt*u*u^T), including its nullspace.
        let change=det::expm1(-norm2/f64::from(b.rate))*along;
        let loss=damp_board(b.board_half_damping.as_ref().unwrap(),&v,&mut out);
        for i in 0..r {assert!((out[i]-v[i]-change*axis[i]).abs()<2e-13);}
        assert!(loss>0.0);
        assert!((loss+0.5*out.iter().map(|x|x*x).sum::<f64>()
            -0.5*v.iter().map(|x|x*x).sum::<f64>()).abs()<1e-15);
        let n=b.modes.len();b.v[n..].copy_from_slice(&v);
        let force=vec![0.05;b.contact_strings.len()];
        for _ in 0..24 {
            let before=b.energy();b.predict();b.finish(&force);
            let work:f64=force.iter().enumerate().map(|(i,f)|
                f*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum();
            assert!((b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work).abs()<1e-11);
            assert!(b.last_modal_loss_j>=-1e-14);b.commit();
        }
    }

    #[test]
    fn dense_board_damping_admission_is_atomic_cold_and_respects_global_lossless() {
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let board=super::super::board::demonstration();let r=board.len();
        for enabled in [false,true] {
            let mut b=Bank::new(&[scale[48]],&board,192_000,21_600.0,4,enabled).unwrap();
            let original=b.board_half_damping.clone();
            let mut asymmetric=identity(r);asymmetric[1]=1.0;
            let mut indefinite=identity(r);indefinite[1]=2.0;indefinite[r]=2.0;
            for bad in [vec![],vec![f64::NAN;r*r],asymmetric,indefinite] {
                assert!(b.configure_bare_board_damping(&bad).is_err());
                assert_eq!(b.board_half_damping,original);assert!(!b.bare_board_damping_configured);
            }
            b.q[0]=1e-6;assert!(b.configure_bare_board_damping(&identity(r)).is_err());b.q[0]=0.0;
            assert_eq!(b.board_half_damping,original);
            let c:Vec<f64>=(0..r*r).map(|k|((k/r+1)*(k%r+1)) as f64).collect();
            let loaded=b.configure_bare_board_damping(&c).unwrap();
            assert_eq!(b.board_half_damping.is_some(),enabled);
            if !enabled {assert!(loaded.iter().all(|x|*x==0.0));}
            let accepted=b.board_half_damping.clone();
            assert!(b.configure_bare_board_damping(&vec![0.0;r*r]).is_err());
            assert_eq!(b.board_half_damping,accepted);
        }
    }

    #[test]
    fn broader_board_retains_high_modes_and_closes_force_work() {
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let board:Vec<BoardMode>=(0..48).map(|i|BoardMode {
            frequency_hz:70.0+25.0*i as f64,damping_ratio:0.002+0.0002*i as f64,
            bridge:[0.01;88],volume:0.03,
        }).collect();
        let mut b=Bank::new(&[scale[48]],&board,192_000,21_600.0,24,true).unwrap();
        assert_eq!(b.board_count,48);
        let force=vec![1.0;b.contact_strings.len()];
        for _ in 0..32 {
            let before=b.energy();b.predict();b.finish(&force);
            let work: f64=(0..force.len()).map(|c|force[c]*(b.contact_position(c,&b.next_q)
                -b.contact_position(c,&b.q))).sum();
            let defect=b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work;
            assert!(defect.abs()<1e-10,"{defect:e}");b.commit();
        }
        assert!(b.v[b.modes.len()+32..].iter().any(|v|v.abs()>1e-12));
    }
    #[test]
    fn higher_bandwidth_boards_preserve_reciprocal_work_above_the_old_limit() {
        let c=super::super::geometry::demonstration_scale().unwrap()[48];
        for count in [33,64,MAX_BOARD_MODES] {
            let board:Vec<_>=(0..count).map(|i| {
                let mut bridge=[0.0;88];bridge[48]=0.04*det::cos(i as f64*0.6)/(count as f64).sqrt();
                BoardMode{frequency_hz:90.0+10.0*i as f64,damping_ratio:0.01,bridge,volume:0.02}
            }).collect();
            let mut b=Bank::new(&[c],&board,192_000,21_600.0,12,true).unwrap();
            assert_eq!(b.board_count,count);
            for k in 0..b.q.len() {b.q[k]=1e-8*det::sin(k as f64);b.v[k]=1e-4*det::cos(k as f64);}
            let forces=vec![0.2;b.contact_strings.len()];
            for _ in 0..20 {
                let before=b.energy();b.predict();b.finish(&forces);
                let work=(0..forces.len()).map(|c|forces[c]*(b.contact_position(c,&b.next_q)-b.contact_position(c,&b.q))).sum::<f64>();
                let defect=b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work;
                assert!(defect.abs()<1e-10,"{count} modes: {defect:e}");
                b.commit();
            }
            let too_many=vec![board[0].clone();MAX_BOARD_MODES+1];
            assert!(Bank::new(&[c],&too_many,192_000,21_600.0,12,true).is_err());
        }
    }

    #[test]
    fn bass_retention_can_reach_the_audio_band_instead_of_stopping_at_128_partials() {
        let mut c=super::super::geometry::demonstration_scale().unwrap()[0];
        c.unison=1;c.duplex_length_m=0.0;
        let mut b=Bank::new(&[c],&super::super::board::demonstration(),192_000,
            21_600.0,MAX_STRING_MODES,false).unwrap();
        assert!(b.modes.len()>128);
        assert!(b.modes.len()<=MAX_STRING_MODES);
        assert!(b.modes.iter().all(|m|m.omega/TAU<=21_600.0));
        assert_eq!(b.capped_string_segments,0);
        let preview=Bank::new(&[c],&super::super::board::demonstration(),192_000,
            21_600.0,24,false).unwrap();
        assert_eq!(preview.modes.len(),24);
        assert_eq!(preview.capped_string_segments,1);
        let forces=[0.1];
        for _ in 0..20 {
            let before=b.energy();b.predict();b.finish(&forces);
            let work=forces[0]*(b.contact_position(0,&b.next_q)-b.contact_position(0,&b.q));
            let defect=b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work;
            assert!(defect.abs()<1e-10,"wide bass work defect {defect:e}");b.commit();
        }
    }

    #[test]
    fn cached_force_response_matches_the_original_solve_for_sparse_and_zero_forces() {
        let mut b=bank(true);let mut force=vec![0.0;b.contact_strings.len()];
        for k in 0..b.q.len(){b.q[k]=1e-8*det::sin(k as f64);b.v[k]=1e-4*det::cos(k as f64);}
        let pointers=(b.next_q.as_ptr(),b.next_v.as_ptr(),b.contact_response.as_ptr());
        for sample in 0..100 {
            force.fill(0.0);if sample%3!=0 {let c=sample%force.len();force[c]=0.2;}
            let (free,q,v,loss)=unfactored(&b,&force);
            b.predict();close(&free,&b.free_contact);b.finish(&force);
            close(&q,&b.next_q);close(&v,&b.next_v);
            assert!((loss-b.last_modal_loss_j).abs()<1e-12);b.commit();
        }
        assert_eq!(pointers,(b.next_q.as_ptr(),b.next_v.as_ptr(),b.contact_response.as_ptr()));
    }
}

#[cfg(test)]
#[path="polarization_bank_tests.rs"]
mod polarization_tests;
