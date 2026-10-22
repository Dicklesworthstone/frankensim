//! The mechanical state a geometry describes (IRP plan §1.4). One owner for
//! imported geometry (fs-io representation admission) and generated parts
//! (fs-mesh solid parts); a consumer that needs it never guesses it.

/// The mechanical state the geometry describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceState {
    /// No preload (stress-free).
    Unloaded,
    /// Under string load.
    Strung,
    /// Assembled, without string load.
    Assembled,
    /// As played.
    Played,
}

impl ReferenceState {
    /// Stable receipt label.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Unloaded => "unloaded",
            Self::Strung => "strung",
            Self::Assembled => "assembled",
            Self::Played => "played",
        }
    }
}
