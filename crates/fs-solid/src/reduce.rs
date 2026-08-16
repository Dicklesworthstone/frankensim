//! Offline lip/reed reduce lab (music bead
//! `frankensim-music-v8-root-3ez8g.3.3`): a 3-D leaflet mesh + a
//! material card become a VALVE CARD — reduced masses, stiffnesses,
//! damping, and rest-orifice geometry measured from the mesh — with the
//! provenance chain and the retained-compliance disclosure that make
//! "that reed" different from "a plausible reed". The mesh itself never
//! runs at audio rate; the card is what
//! `fs_phs::mass_spring_damper` + `bernoulli_volume_flow` islands (and
//! `fs_scenario::BeatingReed`) consume.
//!
//! REDUCTION RECIPE (no new FEM — D23): [`crate::TetLinearElasticProblem`]
//! assembles the ordinary `(K, M)` pencil; [`fs_modal::slice_window`]
//! certifies the low modes; each retained mode is RE-NORMALIZED from the
//! solver's `phi^T M phi = 1` convention to UNIT MEAN TRANSLATION of the
//! orifice face along the opening axis (mass normalization would destroy
//! the units — the recorded fs-solid lesson), giving physical
//! `m_eff = phi^T M phi` [kg] and `k_eff = lambda m_eff` [N/m].
//!
//! RETAINED-COMPLIANCE HONESTY: the modal compliance
//! `sum Gamma_n^2 / lambda_n` of the retained modes is compared against
//! the EXACT static compliance `r^T K^{-1} r` (sparse LDLT solve) for
//! the face-opening load pattern; a reduction whose retained fraction
//! falls below the authored floor REFUSES rather than shipping a stiff
//! lie.
//!
//! DAMPING HONESTY: tissue loss enters as a CARD FIELD (a loss factor,
//! measured when a tissue pack exists, authored + Estimate-labeled until
//! then; `c = eta sqrt(k m)` at resonance). The lab does not invent Q
//! values.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::state_point::IsotropicElasticStatePoint;
use fs_modal::{SliceOptions, slice_window};

use crate::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
};

/// Domain-separated identity for valve cards.
pub const VALVE_CARD_HASH_DOMAIN: &str = "org.frankensim.fs-solid.valve-card.v1";
/// Card schema version (bumps refuse old bytes, never reinterpret).
pub const VALVE_CARD_SCHEMA_VERSION: u32 = 1;

/// Typed refusals from the reduce lab.
#[derive(Debug)]
pub enum ReduceError {
    /// A request parameter is unusable.
    Invalid {
        /// Diagnosis.
        what: &'static str,
    },
    /// The elasticity assembly refused.
    Assembly(crate::linear3::TetElasticError),
    /// The modal solve refused.
    Modal(fs_modal::ModalError),
    /// The static-compliance factorization refused.
    Factor(fs_sparse::LdltError),
    /// The retained modes carry too little of the opening compliance.
    RetainedComplianceTooLow {
        /// Measured retained fraction.
        fraction: f64,
        /// The authored floor.
        floor: f64,
    },
    /// Card bytes failed decode/verify.
    Card {
        /// Diagnosis.
        what: &'static str,
    },
}

impl core::fmt::Display for ReduceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ReduceError::Invalid { what } => write!(f, "FS-SOLID-REDUCE: {what}"),
            ReduceError::Assembly(e) => write!(f, "FS-SOLID-REDUCE assembly: {e:?}"),
            ReduceError::Modal(e) => write!(f, "FS-SOLID-REDUCE modal: {e:?}"),
            ReduceError::Factor(e) => write!(f, "FS-SOLID-REDUCE factor: {e:?}"),
            ReduceError::RetainedComplianceTooLow { fraction, floor } => write!(
                f,
                "FS-SOLID-REDUCE: retained modes carry {fraction:.4} of the opening \
                 compliance, below the authored floor {floor:.4} — a reduction that \
                 discards the compliance refuses rather than shipping a stiff lie"
            ),
            ReduceError::Card { what } => write!(f, "FS-SOLID-REDUCE card: {what}"),
        }
    }
}

impl core::error::Error for ReduceError {}

impl From<crate::linear3::TetElasticError> for ReduceError {
    fn from(e: crate::linear3::TetElasticError) -> Self {
        ReduceError::Assembly(e)
    }
}

impl From<fs_modal::ModalError> for ReduceError {
    fn from(e: fs_modal::ModalError) -> Self {
        ReduceError::Modal(e)
    }
}

impl From<fs_sparse::LdltError> for ReduceError {
    fn from(e: fs_sparse::LdltError) -> Self {
        ReduceError::Factor(e)
    }
}

/// The orifice: which mesh nodes form the moving valve face, which axis
/// opens, and where the opposing surface sits.
#[derive(Debug, Clone)]
pub struct OrificeSpec<'a> {
    /// Node indices of the moving face (the free lip/reed edge surface).
    pub face_nodes: &'a [usize],
    /// Unit opening axis (the face moves along this to close the gap).
    pub opening_axis: [f64; 3],
    /// Unit width axis (the slit's long direction).
    pub width_axis: [f64; 3],
    /// Offset of the opposing surface along the opening axis [m]; the
    /// rest gap is measured FROM THE MESH as the CLOSEST approach
    /// `min(face . axis) - offset`. Face nodes are the surface strip
    /// FACING the opposing plane (the underside of a lip/reed near its
    /// free edge); the third axis (opening x width) is the flow
    /// direction, whose face extent is the channel's effective
    /// thickness.
    pub opposing_plane_offset_m: f64,
}

/// A reduce-lab request.
pub struct ValveCardRequest<'a> {
    /// Leaflet vertex coordinates [m].
    pub nodes_m: &'a [[f64; 3]],
    /// Conforming P1 tetrahedra.
    pub tetrahedra: &'a [[usize; 4]],
    /// Strongly fixed DOFs (`3*node + component`) — the clamped root.
    pub fixed_dofs: &'a [usize],
    /// Resolved material state (E, nu, rho + receipt identity).
    pub material: &'a IsotropicElasticStatePoint,
    /// The orifice geometry spec.
    pub orifice: OrificeSpec<'a>,
    /// Modal window [Hz] handed to `slice_window`.
    pub window_hz: (f64, f64),
    /// Refuse when the retained modal compliance falls below this
    /// fraction of the exact static compliance.
    pub retained_compliance_floor: f64,
    /// AUTHORED structural loss factor eta (Estimate until a measured
    /// tissue/cane pack exists); `c = eta sqrt(k m)`.
    pub loss_factor: f64,
    /// Caller-supplied source identity (mesh + geometry provenance).
    pub source_id: &'a str,
}

/// One reduced mode of the card.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedMode {
    /// Natural frequency [Hz].
    pub frequency_hz: f64,
    /// Effective mass with the face-unit-translation normalization [kg].
    pub mass_kg: f64,
    /// `lambda * mass` [N/m].
    pub stiffness_n_m: f64,
    /// `eta * sqrt(k m)` [N s/m] from the card's authored loss factor.
    pub damping_n_s_m: f64,
    /// Opening-force projection `Gamma = phi^T r` for the unit face
    /// force pattern (phi M-normalized; N/kg^(1/2) per unit pattern).
    pub participation_sqrt_kg: f64,
}

/// The minted valve card.
#[derive(Debug, Clone, PartialEq)]
pub struct ValveCard {
    /// Content identity over the canonical bytes (domain-separated).
    pub identity: ContentHash,
    /// Caller-supplied source identity.
    pub source_id: String,
    /// Material state receipt identity (fs-material resolution chain).
    pub material_identity: ContentHash,
    /// Mesh digest (nodes + tets, bit-exact).
    pub mesh_digest: ContentHash,
    /// Reduced modes, ascending frequency (index 0 drives the 1-DOF
    /// island; 0..2 drive the two-mass realization).
    pub modes: Vec<ReducedMode>,
    /// Rest gap measured from the mesh [m].
    pub rest_gap_m: f64,
    /// Slit width measured from the mesh [m].
    pub width_m: f64,
    /// Effective thickness measured from the mesh [m].
    pub effective_thickness_m: f64,
    /// Retained modal compliance / exact static compliance.
    pub retained_compliance_fraction: f64,
    /// The authored loss factor (Estimate until measured).
    pub loss_factor: f64,
    /// Total leaflet mass [kg].
    pub total_mass_kg: f64,
}

fn unit(v: [f64; 3]) -> Result<[f64; 3], ReduceError> {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if !(n > 0.0 && n.is_finite()) {
        return Err(ReduceError::Invalid {
            what: "axes must be nonzero finite vectors",
        });
    }
    Ok([v[0] / n, v[1] / n, v[2] / n])
}

fn mesh_digest(nodes: &[[f64; 3]], tets: &[[usize; 4]]) -> ContentHash {
    let mut h = DomainHasher::new("org.frankensim.fs-solid.valve-mesh.v1");
    h.update(&(nodes.len() as u64).to_le_bytes());
    for n in nodes {
        for c in n {
            h.update(&c.to_bits().to_le_bytes());
        }
    }
    h.update(&(tets.len() as u64).to_le_bytes());
    for t in tets {
        for &i in t {
            h.update(&(i as u64).to_le_bytes());
        }
    }
    h.finalize()
}

/// Run the reduction: assemble, certify modes, renormalize, disclose the
/// retained compliance, measure the orifice, and mint the card.
///
/// # Errors
/// [`ReduceError`] on any refusal (see the module doc for the
/// retained-compliance law).
#[allow(clippy::too_many_lines)] // one reduction pipeline, kept whole
pub fn reduce_valve(request: &ValveCardRequest<'_>, cx: &Cx<'_>) -> Result<ValveCard, ReduceError> {
    if request.face_nodes_invalid() {
        return Err(ReduceError::Invalid {
            what: "orifice face nodes must be non-empty and in range",
        });
    }
    if !(request.window_hz.0 > 0.0 && request.window_hz.1 > request.window_hz.0) {
        return Err(ReduceError::Invalid {
            what: "modal window needs 0 < lo < hi",
        });
    }
    if !(request.retained_compliance_floor > 0.0 && request.retained_compliance_floor <= 1.0) {
        return Err(ReduceError::Invalid {
            what: "retained-compliance floor must be in (0, 1]",
        });
    }
    if !(request.loss_factor >= 0.0 && request.loss_factor.is_finite()) {
        return Err(ReduceError::Invalid {
            what: "loss factor must be finite and non-negative",
        });
    }
    let axis = unit(request.orifice.opening_axis)?;
    let width_axis = unit(request.orifice.width_axis)?;

    // 1. Assemble the (K, M) pencil.
    let material = TetElasticMaterial::from_resolved_elastic_state(request.material);
    let problem = TetLinearElasticProblem {
        nodes_m: request.nodes_m,
        tetrahedra: request.tetrahedra,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: request.fixed_dofs,
        budget: TetAssemblyBudget::default(),
    };
    let assembly = problem.assemble(cx)?;

    // 2. Certified low modes.
    let tau = core::f64::consts::TAU;
    let window = (
        (tau * request.window_hz.0).powi(2),
        (tau * request.window_hz.1).powi(2),
    );
    let report = slice_window(
        &assembly.stiffness,
        &assembly.mass,
        window,
        &SliceOptions::default(),
    )?;
    if report.modes.is_empty() {
        return Err(ReduceError::Invalid {
            what: "no modes in the requested window",
        });
    }

    // 3. Opening-direction load pattern on the reduced DOFs: the unit
    //    translation of the FACE nodes along the opening axis.
    let n_red = assembly.free_dofs.len();
    let mut r = vec![0.0f64; n_red];
    let mut face_in_free = 0usize;
    for (red, &full) in assembly.free_dofs.iter().enumerate() {
        let node = full / 3;
        let comp = full % 3;
        if request.orifice.face_nodes.contains(&node) {
            r[red] = axis[comp];
            face_in_free += 1;
        }
    }
    if face_in_free == 0 {
        return Err(ReduceError::Invalid {
            what: "every orifice face node is fixed; the valve cannot move",
        });
    }

    // 4. Exact static compliance r^T K^{-1} r (sparse LDLT).
    let symbolic =
        fs_sparse::SymbolicLdlt::analyze(&assembly.stiffness, fs_sparse::DirectOrdering::Amd)?;
    let factor = symbolic.factor(&assembly.stiffness, &fs_sparse::LdltOptions::default())?;
    let k_inv_r = factor.solve(&r);
    let static_compliance: f64 = r.iter().zip(&k_inv_r).map(|(a, b)| a * b).sum();
    if !(static_compliance > 0.0 && static_compliance.is_finite()) {
        return Err(ReduceError::Invalid {
            what: "static opening compliance is not positive finite",
        });
    }

    // 5. Participations + retained-compliance fraction. For M-orthonormal
    //    modes, K^{-1} = sum phi phi^T / lambda, so the modal compliance
    //    of the FORCE pattern r is sum (phi^T r)^2 / lambda — the
    //    projection is phi^T r directly (phi^T M r is the
    //    base-excitation participation, a different animal; the first
    //    build used it and retained 1e-11 of the compliance).
    let mut modes = Vec::new();
    let mut retained = 0.0f64;
    let mean_face = |phi: &[f64]| -> f64 {
        // Mean opening-axis translation of the face nodes.
        let mut acc = 0.0f64;
        let mut count = 0usize;
        for (red, &full) in assembly.free_dofs.iter().enumerate() {
            let node = full / 3;
            let comp = full % 3;
            if request.orifice.face_nodes.contains(&node) {
                acc += phi[red] * axis[comp];
                if comp == 0 {
                    count += 1;
                }
            }
        }
        acc / (count.max(1) as f64)
    };
    for pair in &report.modes {
        let gamma: f64 = pair.phi.iter().zip(&r).map(|(a, b)| a * b).sum();
        retained += gamma * gamma / pair.lambda;
        let u_face = mean_face(&pair.phi);
        if u_face.abs() < 1e-30 {
            // A mode that does not move the face contributes no valve
            // dynamics; keep it out of the card's realizations.
            continue;
        }
        let scale = 1.0 / u_face;
        // phi_hat = phi * scale; m_eff = phi_hat^T M phi_hat = scale^2
        // (solver convention phi^T M phi = 1).
        let mass_kg = scale * scale;
        let stiffness_n_m = pair.lambda * mass_kg;
        let damping_n_s_m = request.loss_factor * (stiffness_n_m * mass_kg).sqrt();
        modes.push(ReducedMode {
            frequency_hz: pair.lambda.sqrt() / tau,
            mass_kg,
            stiffness_n_m,
            damping_n_s_m,
            participation_sqrt_kg: gamma,
        });
    }
    let fraction = retained / static_compliance;
    if fraction < request.retained_compliance_floor {
        return Err(ReduceError::RetainedComplianceTooLow {
            fraction,
            floor: request.retained_compliance_floor,
        });
    }
    if modes.is_empty() {
        return Err(ReduceError::Invalid {
            what: "no retained mode moves the orifice face",
        });
    }

    // 6. Orifice geometry MEASURED from the mesh.
    let thickness_axis = [
        axis[1] * width_axis[2] - axis[2] * width_axis[1],
        axis[2] * width_axis[0] - axis[0] * width_axis[2],
        axis[0] * width_axis[1] - axis[1] * width_axis[0],
    ];
    let mut gap_min = f64::INFINITY;
    let mut w_lo = f64::INFINITY;
    let mut w_hi = f64::NEG_INFINITY;
    let mut t_lo = f64::INFINITY;
    let mut t_hi = f64::NEG_INFINITY;
    for &node in request.orifice.face_nodes {
        let p = request.nodes_m[node];
        let a = p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2];
        let w = p[0] * width_axis[0] + p[1] * width_axis[1] + p[2] * width_axis[2];
        let t = p[0] * thickness_axis[0] + p[1] * thickness_axis[1] + p[2] * thickness_axis[2];
        gap_min = gap_min.min(a - request.orifice.opposing_plane_offset_m);
        w_lo = w_lo.min(w);
        w_hi = w_hi.max(w);
        t_lo = t_lo.min(t);
        t_hi = t_hi.max(t);
    }
    let rest_gap_m = gap_min;
    if !(rest_gap_m > 0.0) {
        return Err(ReduceError::Invalid {
            what: "measured rest gap is not positive (face behind the opposing plane)",
        });
    }
    let width_m = w_hi - w_lo;
    let effective_thickness_m = t_hi - t_lo;

    let mut card = ValveCard {
        identity: ContentHash([0u8; 32]),
        source_id: request.source_id.to_string(),
        material_identity: request.material.resolved().identity(),
        mesh_digest: mesh_digest(request.nodes_m, request.tetrahedra),
        modes,
        rest_gap_m,
        width_m,
        effective_thickness_m,
        retained_compliance_fraction: fraction,
        loss_factor: request.loss_factor,
        total_mass_kg: assembly.total_mass_kg,
    };
    card.identity = card.recomputed_identity();
    Ok(card)
}

impl ValveCardRequest<'_> {
    fn face_nodes_invalid(&self) -> bool {
        self.orifice.face_nodes.is_empty()
            || self
                .orifice
                .face_nodes
                .iter()
                .any(|&n| n >= self.nodes_m.len())
    }
}

impl ValveCard {
    /// The card's identity over its canonical fields (schema-versioned,
    /// domain-separated; upstream identities first, then f64 bits, then
    /// length-prefixed mode payload).
    #[must_use]
    pub fn recomputed_identity(&self) -> ContentHash {
        let mut h = DomainHasher::new(VALVE_CARD_HASH_DOMAIN);
        h.update(&VALVE_CARD_SCHEMA_VERSION.to_le_bytes());
        h.update(self.material_identity.as_bytes());
        h.update(self.mesh_digest.as_bytes());
        h.update(self.source_id.as_bytes());
        h.update(&[0u8]);
        for v in [
            self.rest_gap_m,
            self.width_m,
            self.effective_thickness_m,
            self.retained_compliance_fraction,
            self.loss_factor,
            self.total_mass_kg,
        ] {
            h.update(&v.to_bits().to_le_bytes());
        }
        h.update(&(self.modes.len() as u64).to_le_bytes());
        for m in &self.modes {
            for v in [
                m.frequency_hz,
                m.mass_kg,
                m.stiffness_n_m,
                m.damping_n_s_m,
                m.participation_sqrt_kg,
            ] {
                h.update(&v.to_bits().to_le_bytes());
            }
        }
        h.finalize()
    }

    /// Canonical line-oriented bytes (schema line first, fixed field
    /// order, `{:e}` floats, identity LAST so decode can verify).
    #[must_use]
    pub fn to_canonical_bytes(&self) -> Vec<u8> {
        let mut s = String::new();
        s.push_str("frankensim-valve-card-v1\n");
        s.push_str(&format!("source_id\t{}\n", self.source_id));
        s.push_str(&format!(
            "material_identity\t{}\n",
            self.material_identity.to_hex()
        ));
        s.push_str(&format!("mesh_digest\t{}\n", self.mesh_digest.to_hex()));
        for (key, v) in [
            ("rest_gap_m", self.rest_gap_m),
            ("width_m", self.width_m),
            ("effective_thickness_m", self.effective_thickness_m),
            (
                "retained_compliance_fraction",
                self.retained_compliance_fraction,
            ),
            ("loss_factor", self.loss_factor),
            ("total_mass_kg", self.total_mass_kg),
        ] {
            s.push_str(&format!("{key}\t{v:e}\n"));
        }
        s.push_str(&format!("modes\t{}\n", self.modes.len()));
        for m in &self.modes {
            s.push_str(&format!(
                "mode\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n",
                m.frequency_hz,
                m.mass_kg,
                m.stiffness_n_m,
                m.damping_n_s_m,
                m.participation_sqrt_kg
            ));
        }
        s.push_str(&format!("identity\t{}\n", self.identity.to_hex()));
        s.into_bytes()
    }

    /// Decode + VERIFY canonical bytes (identity recomputed and matched;
    /// tampered bytes refuse).
    ///
    /// # Errors
    /// [`ReduceError::Card`] naming the violated rule.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<ValveCard, ReduceError> {
        let text = core::str::from_utf8(bytes).map_err(|_| ReduceError::Card {
            what: "card bytes are not UTF-8",
        })?;
        let lines: Vec<&str> = text.lines().collect();
        let mut cursor = 0usize;
        if lines.first().copied() != Some("frankensim-valve-card-v1") {
            return Err(ReduceError::Card {
                what: "schema line mismatch",
            });
        }
        cursor += 1;
        fn take<'a>(
            lines: &[&'a str],
            cursor: &mut usize,
            key: &str,
        ) -> Result<&'a str, ReduceError> {
            let line = lines.get(*cursor).copied().ok_or(ReduceError::Card {
                what: "truncated card",
            })?;
            *cursor += 1;
            let (k, v) = line.split_once('\t').ok_or(ReduceError::Card {
                what: "malformed card line",
            })?;
            if k != key {
                return Err(ReduceError::Card {
                    what: "card field order violated",
                });
            }
            Ok(v)
        }
        fn take_num(lines: &[&str], cursor: &mut usize, key: &str) -> Result<f64, ReduceError> {
            take(lines, cursor, key)?
                .parse::<f64>()
                .map_err(|_| ReduceError::Card {
                    what: "bad numeric field",
                })
        }
        let source_id = take(&lines, &mut cursor, "source_id")?.to_string();
        let material_identity =
            ContentHash::from_hex(take(&lines, &mut cursor, "material_identity")?).ok_or(
                ReduceError::Card {
                    what: "bad material identity hex",
                },
            )?;
        let mesh_digest = ContentHash::from_hex(take(&lines, &mut cursor, "mesh_digest")?).ok_or(
            ReduceError::Card {
                what: "bad mesh digest hex",
            },
        )?;
        let rest_gap_m = take_num(&lines, &mut cursor, "rest_gap_m")?;
        let width_m = take_num(&lines, &mut cursor, "width_m")?;
        let effective_thickness_m = take_num(&lines, &mut cursor, "effective_thickness_m")?;
        let retained_compliance_fraction =
            take_num(&lines, &mut cursor, "retained_compliance_fraction")?;
        let loss_factor = take_num(&lines, &mut cursor, "loss_factor")?;
        let total_mass_kg = take_num(&lines, &mut cursor, "total_mass_kg")?;
        let n_modes = take(&lines, &mut cursor, "modes")?
            .parse::<usize>()
            .map_err(|_| ReduceError::Card {
                what: "bad mode count",
            })?;
        let mut modes = Vec::with_capacity(n_modes);
        for _ in 0..n_modes {
            let line = lines.get(cursor).copied().ok_or(ReduceError::Card {
                what: "truncated mode row",
            })?;
            cursor += 1;
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() != 6 || cols[0] != "mode" {
                return Err(ReduceError::Card {
                    what: "malformed mode row",
                });
            }
            let mut vals = [0.0f64; 5];
            for (slot, col) in vals.iter_mut().zip(&cols[1..]) {
                *slot = col.parse::<f64>().map_err(|_| ReduceError::Card {
                    what: "bad mode value",
                })?;
            }
            modes.push(ReducedMode {
                frequency_hz: vals[0],
                mass_kg: vals[1],
                stiffness_n_m: vals[2],
                damping_n_s_m: vals[3],
                participation_sqrt_kg: vals[4],
            });
        }
        let identity = ContentHash::from_hex(take(&lines, &mut cursor, "identity")?).ok_or(
            ReduceError::Card {
                what: "bad identity hex",
            },
        )?;
        let card = ValveCard {
            identity,
            source_id,
            material_identity,
            mesh_digest,
            modes,
            rest_gap_m,
            width_m,
            effective_thickness_m,
            retained_compliance_fraction,
            loss_factor,
            total_mass_kg,
        };
        if card.recomputed_identity() != card.identity {
            return Err(ReduceError::Card {
                what: "identity mismatch: card bytes were altered",
            });
        }
        Ok(card)
    }

    /// The card's 1-DOF island as an `fs_scenario::BeatingReed` (massive
    /// branch: mode 0's mass and stiffness; the closing pressure follows
    /// the runtime's own face convention `P_c = k H / (w * 0.025)` so the
    /// massless and massive branches agree — see
    /// fs-couple `reed_structural`).
    #[must_use]
    pub fn beating_reed(
        &self,
        blowing_pressure_pa: f64,
        attack_s: f64,
    ) -> fs_scenario::BeatingReed {
        let m0 = &self.modes[0];
        fs_scenario::BeatingReed {
            rest_opening_m: self.rest_gap_m,
            width_m: self.width_m,
            closing_pressure_pa: m0.stiffness_n_m * self.rest_gap_m / (self.width_m * 0.025),
            blowing_pressure_pa,
            attack_s,
            mass_kg: m0.mass_kg,
            stiffness_n_m: m0.stiffness_n_m,
        }
    }

    /// JSON-lines log rows: the modal ladder with participations and the
    /// retained-compliance disclosure.
    #[must_use]
    pub fn debug_lines(&self) -> Vec<String> {
        let mut out = vec![format!(
            "{{\"suite\":\"fs-solid\",\"case\":\"valve-card\",\"source_id\":\"{}\",\
             \"rest_gap_m\":{:.4e},\"width_m\":{:.4e},\"thickness_m\":{:.4e},\
             \"retained_compliance_fraction\":{:.4},\"loss_factor\":{:.3},\
             \"total_mass_kg\":{:.4e}}}",
            self.source_id,
            self.rest_gap_m,
            self.width_m,
            self.effective_thickness_m,
            self.retained_compliance_fraction,
            self.loss_factor,
            self.total_mass_kg
        )];
        for (i, m) in self.modes.iter().enumerate() {
            out.push(format!(
                "{{\"suite\":\"fs-solid\",\"case\":\"valve-card-mode\",\"index\":{i},\
                 \"f_hz\":{:.3},\"mass_kg\":{:.4e},\"stiffness_n_m\":{:.4e},\
                 \"damping_n_s_m\":{:.4e},\"participation_sqrt_kg\":{:.4e}}}",
                m.frequency_hz,
                m.mass_kg,
                m.stiffness_n_m,
                m.damping_n_s_m,
                m.participation_sqrt_kg
            ));
        }
        out
    }
}

#[cfg(test)]
mod reduce_tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, resolve_isotropic_elastic_state_point,
    };
    use fs_qty::{Dims, Pressure};

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = fs_exec::CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 0,
                    kernel_id: 77,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn verdict(case: &str, pass: bool, detail: &str) {
        println!(
            "{{\"suite\":\"fs-solid\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
            if pass { "pass" } else { "fail" }
        );
        assert!(pass, "case {case}: {detail}");
    }

    /// Synthetic tissue-ish card (AUTHORED, Estimate: no tissue pack
    /// exists in data/matdb — the recorded upgrade path).
    fn tissue_card(e_pa: f64, rho: f64) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", fs_qty::Density::DIMS, rho),
            ("young_modulus", Pressure::DIMS, e_pa),
            ("poisson_ratio", Dims::NONE, 0.30),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    provenance: Provenance {
                        source: format!("authored tissue-like {name} (Estimate)"),
                        license: "CC0-1.0".to_owned(),
                        artifact: None,
                    },
                    observations: Vec::new(),
                })
                .expect("claim");
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "authored-soft-tissue".to_owned(),
                phase: "solid".to_owned(),
                process: "synthetic".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .expect("card")
    }

    /// Freudenthal 6-tet-per-hex slab: `nx x ny x nz` hexes over
    /// `l x w x t` metres. Returns (nodes, tets).
    fn slab_mesh(
        l: f64,
        w: f64,
        t: f64,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
        let node = |i: usize, j: usize, k: usize| -> usize { (i * (ny + 1) + j) * (nz + 1) + k };
        let mut nodes = Vec::new();
        for i in 0..=nx {
            for j in 0..=ny {
                for k in 0..=nz {
                    nodes.push([
                        l * i as f64 / nx as f64,
                        w * j as f64 / ny as f64,
                        t * k as f64 / nz as f64,
                    ]);
                }
            }
        }
        let mut tets = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let v = [
                        node(i, j, k),
                        node(i + 1, j, k),
                        node(i, j + 1, k),
                        node(i + 1, j + 1, k),
                        node(i, j, k + 1),
                        node(i + 1, j, k + 1),
                        node(i, j + 1, k + 1),
                        node(i + 1, j + 1, k + 1),
                    ];
                    // Freudenthal decomposition around the 0-7 diagonal.
                    for tet in [
                        [v[0], v[1], v[3], v[7]],
                        [v[0], v[3], v[2], v[7]],
                        [v[0], v[2], v[6], v[7]],
                        [v[0], v[6], v[4], v[7]],
                        [v[0], v[4], v[5], v[7]],
                        [v[0], v[5], v[1], v[7]],
                    ] {
                        tets.push(tet);
                    }
                }
            }
        }
        (nodes, tets)
    }

    /// Clamp the x = 0 face; the orifice face is the BOTTOM-surface
    /// strip (z = 0) within one element column of the free tip — the
    /// underside that faces the opposing plane.
    fn clamp_and_tip(nodes: &[[f64; 3]], l: f64, strip_m: f64) -> (Vec<usize>, Vec<usize>) {
        let mut fixed = Vec::new();
        let mut tip = Vec::new();
        for (n, p) in nodes.iter().enumerate() {
            if p[0].abs() < 1e-12 {
                fixed.extend_from_slice(&[3 * n, 3 * n + 1, 3 * n + 2]);
            }
            if p[2].abs() < 1e-12 && p[0] > l - strip_m - 1e-12 {
                tip.push(n);
            }
        }
        (fixed, tip)
    }

    struct Fixture {
        nodes: Vec<[f64; 3]>,
        tets: Vec<[usize; 4]>,
        fixed: Vec<usize>,
        tip: Vec<usize>,
        l: f64,
        w: f64,
        t: f64,
    }

    fn cantilever_fixture() -> Fixture {
        let (l, w, t) = (0.010f64, 0.012f64, 0.004f64);
        let (nodes, tets) = slab_mesh(l, w, t, 6, 6, 3);
        let (fixed, tip) = clamp_and_tip(&nodes, l, l / 6.0);
        Fixture {
            nodes,
            tets,
            fixed,
            tip,
            l,
            w,
            t,
        }
    }

    fn resolve_tissue(e_pa: f64, rho: f64) -> fs_material::state_point::IsotropicElasticStatePoint {
        let card = tissue_card(e_pa, rho);
        let point = QueryPoint::new().with("T", 293.15).expect("point");
        resolve_isotropic_elastic_state_point(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("resolve")
    }

    #[test]
    fn rl_001_cantilever_oracle() {
        with_cx(|cx| {
            let fx = cantilever_fixture();
            let e_pa = 1.0e6;
            let rho = 1050.0;
            let state = resolve_tissue(e_pa, rho);
            let request = ValveCardRequest {
                nodes_m: &fx.nodes,
                tetrahedra: &fx.tets,
                fixed_dofs: &fx.fixed,
                material: &state,
                orifice: OrificeSpec {
                    face_nodes: &fx.tip,
                    opening_axis: [0.0, 0.0, 1.0],
                    width_axis: [0.0, 1.0, 0.0],
                    opposing_plane_offset_m: -0.002,
                },
                window_hz: (20.0, 2000.0),
                retained_compliance_floor: 0.5,
                loss_factor: 0.1,
                source_id: "lab/cantilever-slab/24x12x4mm/v1",
            };
            let card = reduce_valve(&request, cx).expect("reduce");
            for line in card.debug_lines() {
                println!("{line}");
            }
            // Euler-Bernoulli fundamental (thin-beam limit; the thick P1
            // slab deviates and the band is authored from measurement).
            let i_area = fx.w * fx.t * fx.t * fx.t / 12.0;
            let a_area = fx.w * fx.t;
            let f_eb = (1.875f64.powi(2) / core::f64::consts::TAU)
                * (e_pa * i_area / (rho * a_area * fx.l.powi(4))).sqrt();
            let f1 = card.modes[0].frequency_hz;
            let f_ratio = f1 / f_eb;
            // Tip-normalized cantilever fundamental: m_eff ~ m_total/4.
            let mass_ratio = card.modes[0].mass_kg / card.total_mass_kg;
            // k = lambda m identity (construction).
            let lam = (core::f64::consts::TAU * f1).powi(2);
            let k_identity = (card.modes[0].stiffness_n_m - lam * card.modes[0].mass_kg).abs()
                / card.modes[0].stiffness_n_m;
            // Static tip stiffness 3EI/L^3 vs the reduced stiffness.
            let k_static = 3.0 * e_pa * i_area / fx.l.powi(3);
            let k_ratio = card.modes[0].stiffness_n_m / k_static;
            // Orifice geometry measured from the mesh.
            // Bottom strip: gap = 0 - (-0.002); width = full w; channel
            // thickness = the one-element strip length along the flow.
            let geom_ok = (card.width_m - fx.w).abs() < 1e-12
                && (card.effective_thickness_m - fx.l / 6.0).abs() < 1e-9
                && (card.rest_gap_m - 0.002).abs() < 1e-9;
            println!(
                "{{\"suite\":\"fs-solid\",\"case\":\"rl-001-numbers\",\"f1_hz\":{f1:.2},\
                 \"f_eb_hz\":{f_eb:.2},\"f_ratio\":{f_ratio:.3},\"mass_ratio\":{mass_ratio:.4},\
                 \"k_ratio_vs_static\":{k_ratio:.3},\
                 \"retained\":{:.4}}}",
                card.retained_compliance_fraction
            );
            // AUTHORED bands from the executed measure run: the thick
            // slab (L/t = 2.5) sits ABOVE the slender Euler-Bernoulli
            // line (plate-strip stiffening E/(1-nu^2) is +4.8% on f, P1
            // discretization adds the rest; measured +11.4%), the
            // BOTTOM-STRIP-normalized effective mass sits above the
            // slender-beam tip-normalized 0.25 (the strip's mean
            // deflection is below the tip's; measured 0.355), the
            // reduced stiffness exceeds the slender 3EI/L^3 consistently
            // (measured 1.814 = f_ratio^2 * mass_ratio/0.25, the
            // self-consistency identity), and the retained low modes
            // carry essentially the whole opening compliance (0.9805).
            let pass = k_identity < 1e-12
                && geom_ok
                && f_ratio > 1.0
                && f_ratio < 1.35
                && mass_ratio > 0.30
                && mass_ratio < 0.42
                && k_ratio > 1.5
                && k_ratio < 2.2
                && card.retained_compliance_fraction > 0.95;
            verdict(
                "rl-001-cantilever-oracle",
                pass,
                &format!(
                    "f1 {f1:.1} Hz vs EB {f_eb:.1} (ratio {f_ratio:.3}); m_eff/m \
                     {mass_ratio:.3}; k vs 3EI/L^3 ratio {k_ratio:.3}; k=lam*m residual \
                     {k_identity:.2e}; geometry {geom_ok}"
                ),
            );
        });
    }

    #[test]
    fn rl_002_card_drives_the_valve_island() {
        with_cx(|cx| {
            // Lip-ish leaflet: the slab with a parabolic thickness taper
            // (an AUTHORED lip-shaped fixture, honestly labeled; a real
            // lip scan slots through the same seam).
            let fx = cantilever_fixture();
            let mut nodes = fx.nodes.clone();
            for p in &mut nodes {
                let s = p[0] / fx.l;
                p[2] *= 1.0 - 0.4 * s * s;
            }
            let state = resolve_tissue(0.8e6, 1050.0);
            let request = ValveCardRequest {
                nodes_m: &nodes,
                tetrahedra: &fx.tets,
                fixed_dofs: &fx.fixed,
                material: &state,
                orifice: OrificeSpec {
                    face_nodes: &fx.tip,
                    opening_axis: [0.0, 0.0, 1.0],
                    width_axis: [0.0, 1.0, 0.0],
                    opposing_plane_offset_m: -0.0015,
                },
                window_hz: (20.0, 2000.0),
                retained_compliance_floor: 0.5,
                loss_factor: 0.15,
                source_id: "lab/lip-tapered-slab/v1",
            };
            let card = reduce_valve(&request, cx).expect("reduce");
            // Smoke loop: the card's 1-DOF island under a steady blowing
            // pressure through the existing fs-phs machinery.
            let m0 = &card.modes[0];
            let msd = fs_phs::mass_spring_damper(m0.mass_kg, m0.stiffness_n_m, m0.damping_n_s_m)
                .expect("msd admits");
            let face_area = card.width_m * card.effective_thickness_m;
            let p_blow = 500.0f64;
            let rho_air = 1.204;
            let dt = 1.0 / 48_000.0;
            let mut x = vec![0.0f64, 0.0];
            let mut worst_defect = 0.0f64;
            let mut opening = card.rest_gap_m;
            let mut flow = 0.0f64;
            for _ in 0..24_000 {
                let force = p_blow * face_area;
                let record = fs_phs::step(&msd, &x, &[force], dt).expect("step");
                worst_defect = worst_defect.max(record.supply_defect().abs());
                x = record.x;
                opening = (card.rest_gap_m + x[0]).max(0.0);
                flow = fs_phs::bernoulli_volume_flow(card.width_m, opening, p_blow, rho_air);
            }
            // Quasi-static equilibrium of the island: q* = F/k.
            let q_star = p_blow * face_area / m0.stiffness_n_m;
            let q_err = (x[0] - q_star).abs() / q_star.abs().max(1e-30);
            let flow_expected = fs_phs::bernoulli_volume_flow(
                card.width_m,
                (card.rest_gap_m + q_star).max(0.0),
                p_blow,
                rho_air,
            );
            let flow_err = (flow - flow_expected).abs() / flow_expected.max(1e-30);
            let pass = q_err < 1e-3 && flow_err < 1e-3 && worst_defect < 1e-8 && opening > 0.0;
            verdict(
                "rl-002-card-drives-the-island",
                pass,
                &format!(
                    "settled q {:.4e} vs q* {q_star:.4e} (rel {q_err:.2e}); flow rel \
                     {flow_err:.2e}; worst supply defect {worst_defect:.2e}; opening {opening:.4e}",
                    x[0]
                ),
            );
        });
    }

    #[test]
    fn rl_003_refusals() {
        with_cx(|cx| {
            let fx = cantilever_fixture();
            let state = resolve_tissue(1.0e6, 1050.0);
            let base = |face: &'static [usize]| ValveCardRequest {
                nodes_m: &fx.nodes,
                tetrahedra: &fx.tets,
                fixed_dofs: &fx.fixed,
                material: &state,
                orifice: OrificeSpec {
                    face_nodes: face,
                    opening_axis: [0.0, 0.0, 1.0],
                    width_axis: [0.0, 1.0, 0.0],
                    opposing_plane_offset_m: -0.002,
                },
                window_hz: (20.0, 2000.0),
                retained_compliance_floor: 0.5,
                loss_factor: 0.1,
                source_id: "lab/refusals/v1",
            };
            let empty = reduce_valve(&base(&[]), cx);
            let oob = reduce_valve(&base(&[usize::MAX]), cx);
            // Retained-compliance floor: a window that EXCLUDES the
            // compliant fundamental keeps only stiff high modes.
            let mut high = base(Box::leak(fx.tip.clone().into_boxed_slice()));
            high.window_hz = (2500.0, 6000.0);
            high.retained_compliance_floor = 0.5;
            let starved = reduce_valve(&high, cx);
            // A face behind the opposing plane has no positive gap.
            let mut behind = base(Box::leak(fx.tip.clone().into_boxed_slice()));
            behind.orifice.opposing_plane_offset_m = 1.0;
            let negative_gap = reduce_valve(&behind, cx);
            let pass = matches!(empty, Err(ReduceError::Invalid { .. }))
                && matches!(oob, Err(ReduceError::Invalid { .. }))
                && matches!(starved, Err(ReduceError::RetainedComplianceTooLow { .. }))
                && matches!(negative_gap, Err(ReduceError::Invalid { .. }));
            verdict(
                "rl-003-refusals",
                pass,
                &format!(
                    "empty {} oob {} starved {} negative-gap {}",
                    empty.is_err(),
                    oob.is_err(),
                    starved.is_err(),
                    negative_gap.is_err()
                ),
            );
        });
    }

    #[test]
    fn rl_004_card_round_trip_and_tamper() {
        with_cx(|cx| {
            let fx = cantilever_fixture();
            let state = resolve_tissue(1.0e6, 1050.0);
            let request = ValveCardRequest {
                nodes_m: &fx.nodes,
                tetrahedra: &fx.tets,
                fixed_dofs: &fx.fixed,
                material: &state,
                orifice: OrificeSpec {
                    face_nodes: &fx.tip,
                    opening_axis: [0.0, 0.0, 1.0],
                    width_axis: [0.0, 1.0, 0.0],
                    opposing_plane_offset_m: -0.002,
                },
                window_hz: (20.0, 2000.0),
                retained_compliance_floor: 0.5,
                loss_factor: 0.1,
                source_id: "lab/round-trip/v1",
            };
            let card = reduce_valve(&request, cx).expect("reduce");
            let bytes = card.to_canonical_bytes();
            let back = ValveCard::from_canonical_bytes(&bytes).expect("round trip");
            let identical = back == card;
            // Tamper with one stiffness digit: the identity must refuse.
            let text = String::from_utf8(bytes).expect("utf8");
            let tampered = text.replacen("mode\t", "mode\t", 1); // no-op guard
            assert_eq!(tampered, text, "guard");
            let mut altered = text.clone();
            let pos = altered.find("rest_gap_m\t").expect("field");
            altered.replace_range(pos + 11..pos + 12, "9");
            let refused = ValveCard::from_canonical_bytes(altered.as_bytes());
            let beating = card.beating_reed(2000.0, 0.01);
            let reed_ok = beating.mass_kg > 0.0
                && beating.stiffness_n_m > 0.0
                && beating.rest_opening_m == card.rest_gap_m;
            let pass = identical && matches!(refused, Err(ReduceError::Card { .. })) && reed_ok;
            verdict(
                "rl-004-round-trip-and-tamper",
                pass,
                &format!(
                    "round-trip identical {identical}; tampered refused {}; BeatingReed \
                     mint ok {reed_ok}",
                    refused.is_err()
                ),
            );
        });
    }
}
