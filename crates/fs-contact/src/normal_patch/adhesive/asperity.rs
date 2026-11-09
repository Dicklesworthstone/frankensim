//! Independent spherical-asperity adhesion at prescribed plane separation.
//!
//! A bin of height z and multiplicity n sees the ORIGINAL JKR approach z-d,
//! where d is the nominal plane separation. Force, real solid contact area,
//! reversible energy and branch-jump release sum with weight n. Every bin
//! retains its own attachment history: a summit that has never touched cannot
//! acquire the tensile branch of a previously attached neighbour.
//!
//! This is the finite-population counterpart of the Fuller--Tabor approximation.
//! See Greenwood (2017), DOI 10.1007/s11249-017-0938-1, sections 1, 1.1 and 3.
//! Heights, curvatures and multiplicities are explicit geometry or quadrature
//! inputs; Ra, optical roughness and material names do not determine them.
//!
//! Asperities are assumed sufficiently separated that their contact patches
//! neither overlap nor interact elastically. The area-fraction limit is only
//! an additional applicability screen, not proof of either assumption. Local
//! parabolic shape, original JKR envelopes and frozen source conditions remain
//! required. There is no spatial patch merging, background compliance, plastic
//! flattening, shear, viscoelasticity or finite-range precontact attraction.
//! In particular, moderate Tabor parameters can substantially change loading
//! hysteresis: retaining the original JKR model is not an accuracy certificate.

use super::{JkrBranch, JkrError, JkrJumpKind, JkrResponse, JkrState, SourceJkrSphereContact};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::{InputAuthority, TriboError, adhesive::AdhesiveContactState};

/// Fixed stack/work limit. Each bin performs bounded original JKR evaluations.
pub const MAX_JKR_ASPERITY_BINS: usize = 256;

/// One explicitly supplied summit, or one positive height/curvature quadrature bin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrAsperity<'a> {
    /// Original material/interface and spherical-curvature owner.
    pub source: &'a SourceJkrSphereContact,
    /// Summit height relative to the declared nominal plane [m].
    pub height_m: f64,
    /// Number of identical independent summits represented by this bin.
    /// Fractional positive weights are allowed for an explicitly supplied
    /// population quadrature; they are not probability-normalized implicitly.
    pub multiplicity: f64,
}

/// Explicit geometry and model-use limits, not propagated error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrAsperityEnvelope {
    pub minimum_separation_m: f64,
    pub maximum_separation_m: f64,
    /// Fixed nominal footprint containing the represented asperities [m2].
    pub nominal_area_m2: f64,
    /// Admitted real-area/nominal-area ratio, strictly between zero and one.
    pub maximum_contact_area_fraction: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrAsperityError {
    Invalid(&'static str),
    /// The original source's refusal is attributed to its population bin.
    Source {
        index: usize,
        error: JkrError,
    },
    /// Geometry, source, bin order, weights or envelope differ from the owner.
    StateMismatch,
    OutsideSeparation,
    ContactAreaFraction {
        actual: f64,
        limit: f64,
    },
    Unrepresentable(&'static str),
    Cancelled,
}
impl core::fmt::Display for JkrAsperityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "independent JKR asperities: {self:?}")
    }
}
impl std::error::Error for JkrAsperityError {}

/// Extensive endpoint observations of the actual retained contact branches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrAsperityResponse {
    pub separation_m: f64,
    /// Signed sum, positive in compression; tension never removes solid area.
    pub normal_force_n: f64,
    /// Sum of the original Hertz contributions, not the positive part of pressure.
    pub repulsive_normal_force_n: f64,
    pub contact_area_m2: f64,
    pub attached_multiplicity: f64,
    /// dF/d(approach) = -dF/d(separation), on the retained smooth branches.
    /// Absent at an asperity spinodal or a jump exactly at the requested endpoint.
    pub tangent_n_per_m: Option<f64>,
    /// Original elastic plus surface energy relative to all summits detached [J].
    pub reversible_energy_j: f64,
    /// Original synthetic authority survives; all other cases remain estimates.
    pub authority: InputAuthority,
}
impl JkrAsperityResponse {
    /// Capacity input for an explicitly selected common junction-shear law.
    /// This transfers real area and repulsive load; it introduces no shear law,
    /// load-sharing rule or contact-pressure distribution of its own.
    pub fn friction_contact_state(&self) -> Result<AdhesiveContactState, TriboError> {
        AdhesiveContactState::new(self.contact_area_m2, self.repulsive_normal_force_n)
    }
}

/// Sealed accepted branch history. Rewriting public observations cannot rebind it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrAsperityState<const N: usize> {
    identity: ContentHash,
    local: [JkrState; N],
    response: JkrAsperityResponse,
}
impl<const N: usize> JkrAsperityState<N> {
    pub const fn response(&self) -> &JkrAsperityResponse {
        &self.response
    }
    pub const fn local_states(&self) -> &[JkrState; N] {
        &self.local
    }
}

/// One bin's exact fixed-separation event; multiple bins can jump together.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrAsperityJump {
    pub kind: JkrJumpKind,
    /// Nominal plane separation at the event [m], even after an overshoot.
    pub separation_m: f64,
    pub multiplicity: f64,
    /// Weighted original jump energy. Its destination is not automatically heat.
    pub released_energy_j: f64,
}

/// Retryable immutable monotone-separation proposal, with no allocation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrAsperityStep<const N: usize> {
    state: JkrAsperityState<N>,
    /// At most one event per bin, in original BIN order, not chronological order.
    pub jumps: [Option<JkrAsperityJump>; N],
    /// Integral sum(n F_i d(delta_i)) = -integral F_total dd [J].
    pub prescribed_work_j: f64,
    /// Sum of original jump releases [J], distinct from reversible storage.
    pub released_energy_j: f64,
}
impl<const N: usize> JkrAsperityStep<N> {
    pub const fn state(&self) -> JkrAsperityState<N> {
        self.state
    }
    pub const fn response(&self) -> &JkrAsperityResponse {
        &self.state.response
    }
}

/// Frozen original sources and explicit finite population. Construction and
/// evaluation use fixed arrays and O(N) bounded source calls, never allocations.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrAsperityContact<'a, const N: usize> {
    asperities: [JkrAsperity<'a>; N],
    envelope: JkrAsperityEnvelope,
    identity: ContentHash,
}
impl<'a, const N: usize> SourceJkrAsperityContact<'a, N> {
    /// Sources may have distinct summit radii, but must retain the same complete
    /// ordered material/interface state. Each source admits its own JKR regime.
    /// The closest admitted separation conservatively bounds all attached areas.
    pub fn new(
        asperities: [JkrAsperity<'a>; N],
        envelope: JkrAsperityEnvelope,
        gate: &CancelGate,
    ) -> Result<Self, JkrAsperityError> {
        poll(gate)?;
        if N == 0 || N > MAX_JKR_ASPERITY_BINS {
            return Err(JkrAsperityError::Invalid("require 1..=256 asperity bins"));
        }
        if !envelope.minimum_separation_m.is_finite()
            || !envelope.maximum_separation_m.is_finite()
            || envelope.minimum_separation_m >= envelope.maximum_separation_m
            || !(envelope.maximum_separation_m - envelope.minimum_separation_m).is_finite()
            || !envelope.nominal_area_m2.is_finite()
            || envelope.nominal_area_m2 <= 0.0
            || !envelope.maximum_contact_area_fraction.is_finite()
            || envelope.maximum_contact_area_fraction <= 0.0
            || envelope.maximum_contact_area_fraction >= 1.0
        {
            return Err(JkrAsperityError::Invalid("population geometry/envelope"));
        }
        let normal_identity = asperities[0].source.normal_state().identity();
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-asperities.v1");
        hash.update(&(N as u64).to_le_bytes());
        let mut maximum_area = Dd::ZERO;
        let mut total_multiplicity = Dd::ZERO;
        for (index, bin) in asperities.iter().enumerate() {
            poll(gate)?;
            if !bin.height_m.is_finite() || !bin.multiplicity.is_finite() || bin.multiplicity <= 0.0
            {
                return Err(JkrAsperityError::Invalid(
                    "finite height and positive multiplicity required",
                ));
            }
            if bin.source.normal_state().identity() != normal_identity {
                return Err(JkrAsperityError::StateMismatch);
            }
            hash.update(bin.source.identity().as_bytes());
            for x in [bin.height_m, bin.multiplicity] {
                hash.update(&x.to_bits().to_le_bytes());
            }
            let approach = finite(
                bin.height_m - envelope.minimum_separation_m,
                "local approach",
            )?;
            finite(bin.height_m - envelope.maximum_separation_m, "local gap")?;
            bin.source
                .validate_indentation(approach)
                .map_err(|error| source_error(index, error))?;
            if approach >= bin.source.minimum_attached_indentation_m() {
                let local = bin
                    .source
                    .evaluate_attached(approach)
                    .map_err(|error| source_error(index, error))?;
                add(&mut maximum_area, bin.multiplicity, local.contact_area_m2)?;
            }
            add(&mut total_multiplicity, 1.0, bin.multiplicity)?;
        }
        finite(total_multiplicity.to_f64(), "total multiplicity")?;
        area_limit(maximum_area.to_f64(), envelope)?;
        for x in [
            envelope.minimum_separation_m,
            envelope.maximum_separation_m,
            envelope.nominal_area_m2,
            envelope.maximum_contact_area_fraction,
        ] {
            hash.update(&x.to_bits().to_le_bytes());
        }
        poll(gate)?;
        Ok(Self {
            asperities,
            envelope,
            identity: hash.finalize(),
        })
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn asperities(&self) -> &[JkrAsperity<'a>; N] {
        &self.asperities
    }
    pub const fn envelope(&self) -> JkrAsperityEnvelope {
        self.envelope
    }

    /// Explicit initially detached surface. Every summit must have nonpositive
    /// local approach. Exactly touching summits attach on the next advance.
    pub fn detached_state(
        &self,
        separation_m: f64,
        gate: &CancelGate,
    ) -> Result<JkrAsperityState<N>, JkrAsperityError> {
        self.initial_state(separation_m, false, gate)
    }
    /// Explicit already-attached initial branch for EVERY bin. This does not
    /// assert how those bonds formed or fabricate the work of their formation.
    pub fn attached_state(
        &self,
        separation_m: f64,
        gate: &CancelGate,
    ) -> Result<JkrAsperityState<N>, JkrAsperityError> {
        self.initial_state(separation_m, true, gate)
    }
    fn initial_state(
        &self,
        separation_m: f64,
        attached: bool,
        gate: &CancelGate,
    ) -> Result<JkrAsperityState<N>, JkrAsperityError> {
        self.admit(separation_m)?;
        let mut local = [self.asperities[0].source.state(JkrBranch::Detached, 0.0); N];
        let mut totals = Totals::default();
        for (index, bin) in self.asperities.iter().enumerate() {
            poll(gate)?;
            let approach = bin.height_m - separation_m;
            let response = if attached {
                bin.source
                    .evaluate_attached(approach)
                    .map_err(|error| source_error(index, error))?
            } else {
                bin.source
                    .detached_state(approach)
                    .map_err(|error| source_error(index, error))?;
                bin.source.detached_response(approach)
            };
            local[index] = response.state();
            totals.add(bin.multiplicity, &response)?;
        }
        let response = totals.response(separation_m, self.envelope)?;
        poll(gate)?;
        Ok(JkrAsperityState {
            identity: self.identity,
            local,
            response,
        })
    }

    /// One monotone nominal-separation segment; physical reversals require
    /// separate calls. Jumps are localized by the original per-asperity owner,
    /// including when a segment crosses several different summit thresholds.
    /// The caller adopts `state()` only after its receiving mechanics accepts
    /// the result. Cancellation or any local refusal publishes no partial state.
    pub fn advance_monotone(
        &self,
        previous: &JkrAsperityState<N>,
        separation_m: f64,
        gate: &CancelGate,
    ) -> Result<JkrAsperityStep<N>, JkrAsperityError> {
        poll(gate)?;
        if previous.identity != self.identity {
            return Err(JkrAsperityError::StateMismatch);
        }
        self.admit(separation_m)?;
        let mut local = previous.local;
        let mut jumps = [None; N];
        let mut totals = Totals::default();
        let mut work = Dd::ZERO;
        let mut release = Dd::ZERO;
        for (index, bin) in self.asperities.iter().enumerate() {
            poll(gate)?;
            let approach = bin.height_m - separation_m;
            let step = bin
                .source
                .advance_monotone(previous.local[index], approach)
                .map_err(|error| source_error(index, error))?;
            local[index] = step.response.state();
            totals.add(bin.multiplicity, &step.response)?;
            add(&mut work, bin.multiplicity, step.prescribed_work_j)?;
            add(&mut release, bin.multiplicity, step.irreversible_work_j)?;
            if let Some(jump) = step.jump {
                jumps[index] = Some(JkrAsperityJump {
                    kind: jump.kind,
                    separation_m: finite(bin.height_m - jump.indentation_m, "event separation")?,
                    multiplicity: bin.multiplicity,
                    released_energy_j: finite(
                        bin.multiplicity * jump.released_energy_j,
                        "event release",
                    )?,
                });
                if jump.indentation_m == approach {
                    totals.tangent = None;
                }
            }
        }
        let response = totals.response(separation_m, self.envelope)?;
        let state = JkrAsperityState {
            identity: self.identity,
            local,
            response,
        };
        poll(gate)?;
        Ok(JkrAsperityStep {
            state,
            jumps,
            prescribed_work_j: finite(work.to_f64(), "population work")?,
            released_energy_j: finite(release.to_f64(), "population release")?,
        })
    }
    fn admit(&self, separation: f64) -> Result<(), JkrAsperityError> {
        if !separation.is_finite()
            || separation < self.envelope.minimum_separation_m
            || separation > self.envelope.maximum_separation_m
        {
            Err(JkrAsperityError::OutsideSeparation)
        } else {
            Ok(())
        }
    }
}

struct Totals {
    force: Dd,
    repulsive: Dd,
    area: Dd,
    attached: Dd,
    energy: Dd,
    tangent: Option<Dd>,
    synthetic: bool,
}
impl Default for Totals {
    fn default() -> Self {
        Self {
            force: Dd::ZERO,
            repulsive: Dd::ZERO,
            area: Dd::ZERO,
            attached: Dd::ZERO,
            energy: Dd::ZERO,
            tangent: Some(Dd::ZERO),
            synthetic: true,
        }
    }
}
impl Totals {
    fn add(&mut self, weight: f64, response: &JkrResponse) -> Result<(), JkrAsperityError> {
        add(&mut self.force, weight, response.normal_force_n)?;
        add(
            &mut self.repulsive,
            weight,
            response.repulsive_normal_force_n,
        )?;
        add(&mut self.area, weight, response.contact_area_m2)?;
        add(&mut self.energy, weight, response.reversible_energy_j)?;
        if response.branch == JkrBranch::Attached {
            add(&mut self.attached, 1.0, weight)?;
        }
        if let (Some(sum), Some(k)) = (&mut self.tangent, response.tangent_n_per_m) {
            add(sum, weight, k)?;
        } else {
            self.tangent = None;
        }
        self.synthetic &= response.authority == InputAuthority::SyntheticFixture;
        Ok(())
    }
    fn response(
        self,
        separation_m: f64,
        envelope: JkrAsperityEnvelope,
    ) -> Result<JkrAsperityResponse, JkrAsperityError> {
        let area = finite(self.area.to_f64(), "population area")?;
        area_limit(area, envelope)?;
        Ok(JkrAsperityResponse {
            separation_m,
            normal_force_n: finite(self.force.to_f64(), "population force")?,
            repulsive_normal_force_n: finite(
                self.repulsive.to_f64(),
                "population repulsive force",
            )?,
            contact_area_m2: area,
            attached_multiplicity: finite(self.attached.to_f64(), "attached multiplicity")?,
            tangent_n_per_m: self
                .tangent
                .map(|k| finite(k.to_f64(), "population tangent"))
                .transpose()?,
            reversible_energy_j: finite(self.energy.to_f64(), "population energy")?,
            authority: if self.synthetic {
                InputAuthority::SyntheticFixture
            } else {
                InputAuthority::Estimated
            },
        })
    }
}
fn source_error(index: usize, error: JkrError) -> JkrAsperityError {
    JkrAsperityError::Source { index, error }
}
fn finite(x: f64, name: &'static str) -> Result<f64, JkrAsperityError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(JkrAsperityError::Unrepresentable(name))
    }
}
fn add(sum: &mut Dd, weight: f64, value: f64) -> Result<(), JkrAsperityError> {
    *sum = *sum + Dd::from_f64(finite(weight * value, "weighted asperity response")?);
    finite(sum.hi, "population sum")?;
    finite(sum.lo, "population residual")?;
    Ok(())
}
fn area_limit(area: f64, envelope: JkrAsperityEnvelope) -> Result<(), JkrAsperityError> {
    let actual = finite(area / envelope.nominal_area_m2, "contact area fraction")?;
    if actual > envelope.maximum_contact_area_fraction {
        Err(JkrAsperityError::ContactAreaFraction {
            actual,
            limit: envelope.maximum_contact_area_fraction,
        })
    } else {
        Ok(())
    }
}
fn poll(gate: &CancelGate) -> Result<(), JkrAsperityError> {
    if gate.is_requested() {
        Err(JkrAsperityError::Cancelled)
    } else {
        Ok(())
    }
}
