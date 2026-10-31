//! Resolve one physical owner from the existing nested dielectric stack.
//!
//! A bound grid replaces the host's complete homogeneous extinction law. Even
//! an empty grid retains that ownership; normalizing it to an absent binding
//! would accidentally restore absorption from the glass material. Regions are
//! selected without allocating or copying the fixed-capacity boundary stack.

use super::*;

#[derive(Clone, Copy)]
pub(super) struct MediumRegion<'view, 'field> {
    pub(super) volume: Option<PathMedium<'view, 'field>>,
    pub(super) host: Option<DielectricGlass>,
    pub(super) replaces_host_extinction: bool,
}

impl<'view, 'field> MediumRegion<'view, 'field> {
    fn from_entry(
        entry: Option<&MediumEntry<'_>>,
        ambient: Option<PathMedium<'view, 'field>>,
        media: Option<&SceneMedia<'view, 'field>>,
    ) -> Self {
        match entry {
            Some(entry) => {
                let volume = media.and_then(|media| media.interior_path(entry.boundary_primitive));
                Self {
                    volume,
                    host: Some(*entry.glass),
                    replaces_host_extinction: volume.is_some(),
                }
            }
            None => Self {
                volume: ambient,
                host: None,
                replaces_host_extinction: false,
            },
        }
    }

    pub(super) fn incident(
        stack: &MediumStack<'_>,
        ambient: Option<PathMedium<'view, 'field>>,
        media: Option<&SceneMedia<'view, 'field>>,
    ) -> Self {
        Self::from_entry(stack.last(), ambient, media)
    }

    pub(super) fn transmitted(
        stack: &MediumStack<'_>,
        transition: MediumTransition<'_>,
        ambient: Option<PathMedium<'view, 'field>>,
        media: Option<&SceneMedia<'view, 'field>>,
    ) -> Self {
        match transition {
            MediumTransition::Enter(entry) => Self::from_entry(Some(&entry), ambient, media),
            MediumTransition::Exit { .. } => Self::from_entry(
                stack
                    .len()
                    .checked_sub(2)
                    .and_then(|index| stack.get(index)),
                ambient,
                media,
            ),
        }
    }

    /// The residual homogeneous law is absent when an explicit grid owns it.
    pub(super) fn homogeneous_host(self) -> Option<DielectricGlass> {
        if self.replaces_host_extinction {
            None
        } else {
            self.host
        }
    }

    pub(super) fn has_source_proposal(self, strategy: DirectStrategy) -> bool {
        strategy != DirectStrategy::BsdfOnly
            && self
                .volume
                .is_some_and(|volume| volume.emission_grid().has_emission_proposal())
    }
}
