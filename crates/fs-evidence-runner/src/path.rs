//! Exact logical-path and ContentStore-key validation for Runner V2.
//!
//! Values in this module are logical names only.  They never contain or imply
//! a physical root, descriptor, handle, credential, transaction, generation,
//! attempt locator, or acquired capability.

use core::cmp::Ordering;

use crate::catalog::PlatformPathProfileV2;

/// Maximum UTF-8 byte length of a logical path or object key.
pub const LOGICAL_PATH_MAX_BYTES: usize = 240;
/// Maximum segment count of a logical path or object key.
pub const LOGICAL_PATH_MAX_SEGMENTS: usize = 32;
/// Reserved first-segment prefix for FrankenSim ContentStore internals.
pub const CONTENT_STORE_FRANKENSIM_RESERVED_PREFIX: &str = "__frankensim_";
/// Reserved first-segment prefix for Runner ContentStore internals.
pub const CONTENT_STORE_RUNNER_RESERVED_PREFIX: &str = "__runner_";

/// Deterministic single-path construction failures.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PathError {
    /// The complete path was empty.
    Empty,
    /// The complete path began with `/`.
    Absolute,
    /// A segment began with an ASCII drive designator such as `C:`.
    DriveDesignator,
    /// A backslash appeared anywhere in the exact input.
    Backslash {
        /// Zero-based byte offset.
        index: usize,
    },
    /// A NUL byte appeared anywhere in the exact input.
    Nul {
        /// Zero-based byte offset.
        index: usize,
    },
    /// A leading, trailing, or doubled slash exposed an empty segment.
    EmptySegment {
        /// Zero-based segment ordinal.
        segment: usize,
    },
    /// A `.` segment appeared.
    DotSegment {
        /// Zero-based segment ordinal.
        segment: usize,
    },
    /// A `..` segment appeared.
    DotDotSegment {
        /// Zero-based segment ordinal.
        segment: usize,
    },
    /// The complete UTF-8 byte length exceeded the frozen limit.
    TooManyBytes {
        /// Observed byte length.
        observed: usize,
        /// Maximum admitted byte length.
        maximum: usize,
    },
    /// The segment count exceeded the frozen limit.
    TooManySegments {
        /// Observed segment count.
        observed: usize,
        /// Maximum admitted segment count.
        maximum: usize,
    },
    /// A ContentStore key's first segment used a reserved prefix.
    ReservedContentStorePrefix {
        /// The exact reserved prefix that matched.
        prefix: &'static str,
    },
}

/// A validated, exact UTF-8, slash-separated bundle-relative logical path.
///
/// ```
/// use fs_evidence_runner::LogicalBundlePathV1;
///
/// let path = LogicalBundlePathV1::new("artifact/result.bin").unwrap();
/// assert_eq!(path.as_str(), "artifact/result.bin");
/// assert_eq!(path.segment_count(), 2);
/// ```
///
/// A raw string cannot be passed where validation is required:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::LogicalBundlePathV1;
///
/// fn consume_validated(_: LogicalBundlePathV1) {}
///
/// consume_validated("../unvalidated".to_owned());
/// ```
///
/// The tuple field is private, so callers cannot mint the nominal wrapper:
///
/// ```compile_fail,E0423
/// use fs_evidence_runner::LogicalBundlePathV1;
///
/// let _unchecked = LogicalBundlePathV1("../unvalidated".to_owned());
/// ```
///
/// The checked path cannot be mutated after validation:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::LogicalBundlePathV1;
///
/// let mut path = LogicalBundlePathV1::new("artifact/result.bin").unwrap();
/// path.0.push_str("/../escape");
/// ```
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct LogicalBundlePathV1(String);

impl LogicalBundlePathV1 {
    /// Validates `value` without normalization or case folding.
    pub fn new(value: &str) -> Result<Self, PathError> {
        validate_logical_path(value)?;
        Ok(Self(value.to_owned()))
    }

    /// Returns the exact validated UTF-8 input.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Returns the exact validated UTF-8 bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_bytes()
    }

    /// Iterates segments in their exact order without allocating.
    #[must_use]
    pub fn segments(&self) -> impl DoubleEndedIterator<Item = &str> {
        self.0.split('/')
    }

    /// Returns the validated segment count.
    #[must_use]
    pub fn segment_count(&self) -> usize {
        self.segments().count()
    }
}

impl Ord for LogicalBundlePathV1 {
    fn cmp(&self, other: &Self) -> Ordering {
        compare_segment_sequences(self.as_str(), other.as_str())
    }
}

impl PartialOrd for LogicalBundlePathV1 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A validated ContentStore logical object key.
///
/// This is nominally distinct from [`LogicalBundlePathV1`] and additionally
/// rejects both reserved first-segment prefixes.
///
/// ```
/// use fs_evidence_runner::ContentStoreObjectKeyV1;
///
/// let key = ContentStoreObjectKeyV1::new("objects/result").unwrap();
/// assert_eq!(key.as_str(), "objects/result");
/// assert_eq!(key.segment_count(), 2);
/// ```
///
/// A validated bundle path is not silently promoted to a ContentStore key:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{ContentStoreObjectKeyV1, LogicalBundlePathV1};
///
/// let path = LogicalBundlePathV1::new("artifact/object").unwrap();
/// let _key: ContentStoreObjectKeyV1 = path;
/// ```
///
/// A checked object key cannot be widened into a reserved namespace after
/// construction:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::ContentStoreObjectKeyV1;
///
/// let mut key = ContentStoreObjectKeyV1::new("objects/result").unwrap();
/// key.0.insert_str(0, "__runner_");
/// ```
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct ContentStoreObjectKeyV1(String);

impl ContentStoreObjectKeyV1 {
    /// Validates `value` exactly and applies ContentStore reserved-prefix
    /// rules.
    pub fn new(value: &str) -> Result<Self, PathError> {
        validate_logical_path(value)?;
        let first = value.split('/').next().ok_or(PathError::Empty)?;
        if first.starts_with(CONTENT_STORE_FRANKENSIM_RESERVED_PREFIX) {
            return Err(PathError::ReservedContentStorePrefix {
                prefix: CONTENT_STORE_FRANKENSIM_RESERVED_PREFIX,
            });
        }
        if first.starts_with(CONTENT_STORE_RUNNER_RESERVED_PREFIX) {
            return Err(PathError::ReservedContentStorePrefix {
                prefix: CONTENT_STORE_RUNNER_RESERVED_PREFIX,
            });
        }
        Ok(Self(value.to_owned()))
    }

    /// Returns the exact validated UTF-8 input.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Returns the exact validated UTF-8 bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_bytes()
    }

    /// Iterates segments in their exact order without allocating.
    #[must_use]
    pub fn segments(&self) -> impl DoubleEndedIterator<Item = &str> {
        self.0.split('/')
    }

    /// Returns the validated segment count.
    #[must_use]
    pub fn segment_count(&self) -> usize {
        self.segments().count()
    }
}

impl Ord for ContentStoreObjectKeyV1 {
    fn cmp(&self, other: &Self) -> Ordering {
        compare_segment_sequences(self.as_str(), other.as_str())
    }
}

impl PartialOrd for ContentStoreObjectKeyV1 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn validate_logical_path(value: &str) -> Result<(), PathError> {
    if value.is_empty() {
        return Err(PathError::Empty);
    }
    if value.len() > LOGICAL_PATH_MAX_BYTES {
        return Err(PathError::TooManyBytes {
            observed: value.len(),
            maximum: LOGICAL_PATH_MAX_BYTES,
        });
    }
    if value.starts_with('/') {
        return Err(PathError::Absolute);
    }
    if let Some(index) = value.as_bytes().iter().position(|byte| *byte == b'\\') {
        return Err(PathError::Backslash { index });
    }
    if let Some(index) = value.as_bytes().iter().position(|byte| *byte == 0) {
        return Err(PathError::Nul { index });
    }

    // The count is deliberately capped at one past the admitted maximum. That
    // keeps work and arithmetic bounded even if the byte ceiling changes in a
    // later schema generation.
    let segment_count = value.split('/').take(LOGICAL_PATH_MAX_SEGMENTS + 1).count();
    if segment_count > LOGICAL_PATH_MAX_SEGMENTS {
        return Err(PathError::TooManySegments {
            observed: segment_count,
            maximum: LOGICAL_PATH_MAX_SEGMENTS,
        });
    }
    for (segment, component) in value.split('/').enumerate() {
        if is_ascii_drive_designator(component.as_bytes()) {
            return Err(PathError::DriveDesignator);
        }
        match component {
            "" => return Err(PathError::EmptySegment { segment }),
            "." => return Err(PathError::DotSegment { segment }),
            ".." => return Err(PathError::DotDotSegment { segment }),
            _ => {}
        }
    }
    Ok(())
}

fn is_ascii_drive_designator(bytes: &[u8]) -> bool {
    matches!(bytes, [letter, b':', ..] if letter.is_ascii_alphabetic())
}

fn compare_segment_sequences(left: &str, right: &str) -> Ordering {
    let mut left_segments = left.split('/');
    let mut right_segments = right.split('/');
    loop {
        match (left_segments.next(), right_segments.next()) {
            (Some(left_segment), Some(right_segment)) => {
                let ordering = left_segment.as_bytes().cmp(right_segment.as_bytes());
                if ordering != Ordering::Equal {
                    return ordering;
                }
            }
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (None, None) => return Ordering::Equal,
        }
    }
}

fn is_strict_segment_prefix(prefix: &str, descendant: &str) -> bool {
    let mut prefix_segments = prefix.split('/');
    let mut descendant_segments = descendant.split('/');
    let mut matched = 0_usize;

    loop {
        match (prefix_segments.next(), descendant_segments.next()) {
            (Some(left), Some(right)) if left.as_bytes() == right.as_bytes() => {
                matched += 1;
            }
            (Some(_), Some(_) | None) | (None, None) => return false,
            (None, Some(_)) => return matched > 0,
        }
    }
}

fn ascii_folded_relation(left: &str, right: &str) -> FoldedRelation {
    let mut left_segments = left.split('/');
    let mut right_segments = right.split('/');
    loop {
        match (left_segments.next(), right_segments.next()) {
            (Some(left_segment), Some(right_segment)) => {
                if !left_segment.bytes().eq(right_segment.bytes())
                    && !left_segment
                        .bytes()
                        .map(|byte| byte.to_ascii_lowercase())
                        .eq(right_segment.bytes().map(|byte| byte.to_ascii_lowercase()))
                {
                    return FoldedRelation::Distinct;
                }
            }
            (None, Some(_)) => return FoldedRelation::StrictPrefix,
            (Some(_), None) => return FoldedRelation::StrictDescendant,
            (None, None) => return FoldedRelation::Equal,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FoldedRelation {
    Distinct,
    Equal,
    StrictPrefix,
    StrictDescendant,
}

/// Deterministic result of locally adjudicating a validated path set.
///
/// Returned path references always use canonical segment-sequence byte order,
/// so changing caller input order does not change the reported pair.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PathSetAdjudicationV1<'a> {
    /// The profile's collision cell was fully decidable and no collision was
    /// found.
    Exact,
    /// Two entries had exactly identical UTF-8 bytes.
    Duplicate {
        /// The duplicated exact path.
        path: &'a str,
    },
    /// One path was a proper segment prefix of another.
    StrictSegmentPrefix {
        /// The shorter exact path.
        prefix: &'a str,
        /// The longer exact path.
        descendant: &'a str,
    },
    /// Distinct exact ASCII paths alias, or form an aliasing segment-prefix
    /// pair, under Windows ASCII case folding.
    WindowsAsciiAlias {
        /// Canonically earlier exact path.
        first: &'a str,
        /// Canonically later exact path.
        second: &'a str,
    },
    /// At least one Windows path contains non-ASCII bytes, whose alias key is
    /// deliberately outside this base slice.
    UnsupportedWindowsNonAsciiAlias {
        /// Canonically first path requiring platform-owned adjudication.
        path: &'a str,
    },
}

impl PathSetAdjudicationV1<'_> {
    /// Returns true only for a fully decidable, collision-free set.
    #[must_use]
    pub const fn is_exact(self) -> bool {
        matches!(self, Self::Exact)
    }
}

/// Adjudicates a set of validated logical bundle paths for one platform
/// profile.
///
/// Posix and ContentStore cells compare exact bytes.  Windows additionally
/// rejects locally decidable ASCII aliases and reports non-ASCII alias
/// adjudication as unsupported.
#[must_use]
pub fn adjudicate_logical_bundle_path_set(
    profile: PlatformPathProfileV2,
    paths: &[LogicalBundlePathV1],
) -> PathSetAdjudicationV1<'_> {
    let strings = paths.iter().map(LogicalBundlePathV1::as_str);
    adjudicate_strings(profile, &strings)
}

/// Adjudicates ContentStore keys under its exact-byte collision cell.
#[must_use]
pub fn adjudicate_content_store_object_key_set(
    keys: &[ContentStoreObjectKeyV1],
) -> PathSetAdjudicationV1<'_> {
    let strings = keys.iter().map(ContentStoreObjectKeyV1::as_str);
    adjudicate_strings(PlatformPathProfileV2::ContentStoreObjectKeyV1, &strings)
}

fn adjudicate_strings<'a, I>(profile: PlatformPathProfileV2, paths: &I) -> PathSetAdjudicationV1<'a>
where
    I: Clone + Iterator<Item = &'a str>,
{
    let mut duplicate: Option<&'a str> = None;
    let mut prefix_pair: Option<(&'a str, &'a str)> = None;
    let mut windows_alias_pair: Option<(&'a str, &'a str)> = None;
    let mut windows_non_ascii: Option<&'a str> = None;

    for (left_index, left) in (*paths).clone().enumerate() {
        if matches!(profile, PlatformPathProfileV2::WindowsHandleRelativeV1)
            && !left.is_ascii()
            && windows_non_ascii
                .is_none_or(|current| compare_segment_sequences(left, current) == Ordering::Less)
        {
            windows_non_ascii = Some(left);
        }

        for right in (*paths).clone().skip(left_index + 1) {
            let (first, second) = if compare_segment_sequences(left, right) == Ordering::Greater {
                (right, left)
            } else {
                (left, right)
            };

            if left.as_bytes() == right.as_bytes() {
                if duplicate.is_none_or(|current| {
                    compare_segment_sequences(first, current) == Ordering::Less
                }) {
                    duplicate = Some(first);
                }
                continue;
            }

            let exact_prefix_pair = if is_strict_segment_prefix(left, right) {
                Some((left, right))
            } else if is_strict_segment_prefix(right, left) {
                Some((right, left))
            } else {
                None
            };
            if let Some(candidate) = exact_prefix_pair {
                if pair_is_earlier(candidate, prefix_pair) {
                    prefix_pair = Some(candidate);
                }
                continue;
            }

            if matches!(profile, PlatformPathProfileV2::WindowsHandleRelativeV1)
                && left.is_ascii()
                && right.is_ascii()
                && !matches!(ascii_folded_relation(left, right), FoldedRelation::Distinct)
                && pair_is_earlier((first, second), windows_alias_pair)
            {
                windows_alias_pair = Some((first, second));
            }
        }
    }

    if let Some(path) = duplicate {
        return PathSetAdjudicationV1::Duplicate { path };
    }
    if let Some((prefix, descendant)) = prefix_pair {
        return PathSetAdjudicationV1::StrictSegmentPrefix { prefix, descendant };
    }
    if let Some((first, second)) = windows_alias_pair {
        return PathSetAdjudicationV1::WindowsAsciiAlias { first, second };
    }
    if let Some(path) = windows_non_ascii {
        return PathSetAdjudicationV1::UnsupportedWindowsNonAsciiAlias { path };
    }
    PathSetAdjudicationV1::Exact
}

fn pair_is_earlier(candidate: (&str, &str), current: Option<(&str, &str)>) -> bool {
    let Some(current) = current else {
        return true;
    };
    match compare_segment_sequences(candidate.0, current.0) {
        Ordering::Less => true,
        Ordering::Greater => false,
        Ordering::Equal => compare_segment_sequences(candidate.1, current.1) == Ordering::Less,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_utf8_bytes_are_preserved_without_normalization() {
        let composed = LogicalBundlePathV1::new("artifact/café").expect("valid");
        let decomposed = LogicalBundlePathV1::new("artifact/cafe\u{301}").expect("valid");
        assert_eq!(composed.as_str().as_bytes(), "artifact/café".as_bytes());
        assert_eq!(
            decomposed.as_str().as_bytes(),
            "artifact/cafe\u{301}".as_bytes()
        );
        assert_ne!(composed, decomposed);
    }

    #[test]
    fn byte_and_segment_boundaries_are_exact() {
        assert!(LogicalBundlePathV1::new(&"a".repeat(240)).is_ok());
        assert_eq!(
            LogicalBundlePathV1::new(&"a".repeat(241)),
            Err(PathError::TooManyBytes {
                observed: 241,
                maximum: 240
            })
        );

        let thirty_two = vec!["a"; 32].join("/");
        let thirty_three = vec!["a"; 33].join("/");
        assert_eq!(
            LogicalBundlePathV1::new(&thirty_two)
                .expect("valid")
                .segment_count(),
            32
        );
        assert_eq!(
            LogicalBundlePathV1::new(&thirty_three),
            Err(PathError::TooManySegments {
                observed: 33,
                maximum: 32
            })
        );
    }

    #[test]
    fn unsafe_or_ambiguous_single_path_forms_refuse() {
        let cases = [
            ("", PathError::Empty),
            ("/absolute", PathError::Absolute),
            ("C:drive", PathError::DriveDesignator),
            ("z:/drive", PathError::DriveDesignator),
            ("nested/C:drive", PathError::DriveDesignator),
            ("a\\b", PathError::Backslash { index: 1 }),
            ("a\0b", PathError::Nul { index: 1 }),
            ("a//b", PathError::EmptySegment { segment: 1 }),
            ("a/", PathError::EmptySegment { segment: 1 }),
            ("./a", PathError::DotSegment { segment: 0 }),
            ("a/.", PathError::DotSegment { segment: 1 }),
            ("../a", PathError::DotDotSegment { segment: 0 }),
            ("a/..", PathError::DotDotSegment { segment: 1 }),
        ];
        for (input, expected) in cases {
            assert_eq!(LogicalBundlePathV1::new(input), Err(expected), "{input:?}");
        }
    }

    #[test]
    fn content_store_rejects_both_exact_reserved_first_segment_prefixes() {
        for (input, prefix) in [
            ("__frankensim_internal/object", "__frankensim_"),
            ("__runner_state/object", "__runner_"),
        ] {
            assert_eq!(
                ContentStoreObjectKeyV1::new(input),
                Err(PathError::ReservedContentStorePrefix { prefix })
            );
        }
        assert!(ContentStoreObjectKeyV1::new("safe/__runner_state").is_ok());
        assert!(ContentStoreObjectKeyV1::new("__Runner_state/object").is_ok());
    }

    #[test]
    fn canonical_order_is_segment_sequence_byte_order() {
        let mut paths = [
            LogicalBundlePathV1::new("a-").expect("valid"),
            LogicalBundlePathV1::new("a/b").expect("valid"),
            LogicalBundlePathV1::new("a").expect("valid"),
            LogicalBundlePathV1::new("b").expect("valid"),
        ];
        paths.sort();
        assert_eq!(
            paths.map(|path| path.as_str().to_owned()),
            ["a", "a/b", "a-", "b"]
        );
    }

    #[test]
    fn duplicate_and_strict_segment_prefix_are_distinct_and_deterministic() {
        let duplicate = [
            LogicalBundlePathV1::new("z").expect("valid"),
            LogicalBundlePathV1::new("a/b").expect("valid"),
            LogicalBundlePathV1::new("a/b").expect("valid"),
        ];
        assert_eq!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                &duplicate
            ),
            PathSetAdjudicationV1::Duplicate { path: "a/b" }
        );

        let prefix = [
            LogicalBundlePathV1::new("a/bb").expect("valid"),
            LogicalBundlePathV1::new("a/b/c").expect("valid"),
            LogicalBundlePathV1::new("a/b").expect("valid"),
        ];
        assert_eq!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                &prefix
            ),
            PathSetAdjudicationV1::StrictSegmentPrefix {
                prefix: "a/b",
                descendant: "a/b/c"
            }
        );
        assert_eq!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                &prefix[..1]
            ),
            PathSetAdjudicationV1::Exact
        );
    }

    #[test]
    fn posix_and_content_store_cells_are_exact_bytewise() {
        let paths = [
            LogicalBundlePathV1::new("Artifact/value").expect("valid"),
            LogicalBundlePathV1::new("artifact/value").expect("valid"),
            LogicalBundlePathV1::new("artifact/café").expect("valid"),
        ];
        assert!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                &paths
            )
            .is_exact()
        );

        let keys = [
            ContentStoreObjectKeyV1::new("Artifact/value").expect("valid"),
            ContentStoreObjectKeyV1::new("artifact/value").expect("valid"),
            ContentStoreObjectKeyV1::new("artifact/café").expect("valid"),
        ];
        assert!(adjudicate_content_store_object_key_set(&keys).is_exact());
    }

    #[test]
    fn windows_ascii_aliases_refuse_and_non_ascii_is_explicitly_unsupported() {
        let aliases = [
            LogicalBundlePathV1::new("Artifact/Value").expect("valid"),
            LogicalBundlePathV1::new("artifact/value").expect("valid"),
        ];
        assert_eq!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                &aliases
            ),
            PathSetAdjudicationV1::WindowsAsciiAlias {
                first: "Artifact/Value",
                second: "artifact/value"
            }
        );

        let alias_prefix = [
            LogicalBundlePathV1::new("Artifact").expect("valid"),
            LogicalBundlePathV1::new("artifact/value").expect("valid"),
        ];
        assert!(matches!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                &alias_prefix
            ),
            PathSetAdjudicationV1::WindowsAsciiAlias { .. }
        ));

        let unicode = [LogicalBundlePathV1::new("artifact/café").expect("valid")];
        assert_eq!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                &unicode
            ),
            PathSetAdjudicationV1::UnsupportedWindowsNonAsciiAlias {
                path: "artifact/café"
            }
        );
    }

    #[test]
    fn set_adjudication_is_invariant_to_input_permutation() {
        let first = [
            LogicalBundlePathV1::new("z").expect("valid"),
            LogicalBundlePathV1::new("A").expect("valid"),
            LogicalBundlePathV1::new("a").expect("valid"),
        ];
        let second = [first[2].clone(), first[0].clone(), first[1].clone()];
        assert_eq!(
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                &first
            ),
            adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                &second
            )
        );
    }
}
