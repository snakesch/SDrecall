//! Crate-wide error type — the single error enum every SDrecall stage crate returns.
//!
//! Ports the scattered `raise RuntimeError(...)` / `KeyError` / `assert` failure
//! modes of the Python pipeline into one typed enum (the T0 interface contract).
//! Variants are added as stages land; the set below is the shared core. All fields
//! are owned `String`/`usize`/`io::Error`, so the enum is self-contained (no
//! dependency on the geometry / path types) and cheap to construct at the boundary.

use thiserror::Error;

/// The SDrecall error enum. See [`Result`] for the crate-wide alias.
#[derive(Error, Debug)]
pub enum SdError {
    /// Filesystem error, tagged with the path it happened on.
    #[error("I/O error on {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },

    /// Any error surfaced from `rust-htslib` (its `Error` is not cleanly
    /// `Send`/`Clone`, so we stringify it at the boundary).
    #[error("htslib error: {0}")]
    Htslib(String),

    /// A CIGAR contained an `M` op where `=`/`X` were required (align with
    /// `minimap2 --eqx`). Carries the offending read name.
    #[error("CIGAR contains an 'M' op (need '='/'X'; align with minimap2 --eqx) on read {qname}")]
    CigarMOp { qname: String },

    /// A BAM read name was not valid UTF-8.
    #[error("read name is not valid UTF-8")]
    NonUtf8ReadName,

    /// Could not infer the genome assembly from the given SD map + reference.
    #[error("cannot infer assembly from sd_map={sd_map} reference={reference}")]
    UnknownAssembly { sd_map: String, reference: String },

    /// A realignment-group label did not match the expected `RG<n>` shape.
    #[error("invalid realignment-group label: {0}")]
    InvalidRgLabel(String),

    /// A BED line failed to parse.
    #[error("BED parse error at line {line}: {msg}")]
    BedParse { line: usize, msg: String },

    /// Any error surfaced while reading/writing VCF/BCF.
    #[error("VCF error: {0}")]
    Vcf(String),

    /// Any error surfaced while reading/writing GraphML.
    #[error("GraphML error: {0}")]
    GraphMl(String),

    /// A numeric / statistics computation failed (e.g. an out-of-domain
    /// distribution parameter).
    #[error("computation error: {0}")]
    Compute(String),

    /// Not enough usable read pairs to compute a statistic; carries the count seen.
    #[error("insufficient usable reads: {0}")]
    InsufficientPairs(usize),
}

/// Crate-wide result alias: every fallible SDrecall API returns `Result<T>`.
pub type Result<T> = std::result::Result<T, SdError>;
