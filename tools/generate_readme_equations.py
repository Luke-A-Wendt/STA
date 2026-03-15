#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / 'README.txt'
GENERATED_DIR = ROOT / 'manuscript' / 'generated'
OUT = GENERATED_DIR / 'readme_equation_catalog.tex'
FULL_OUT = GENERATED_DIR / 'readme_full_conversion.tex'
APPENDIX_OUT = GENERATED_DIR / 'extended_notes_appendix.tex'
BODY_OUT = GENERATED_DIR / 'integrated_notes_body.tex'
GROUP_OUT_DIR = GENERATED_DIR / 'integrated_sections'

# Unicode sub/superscript maps.
SUB_MAP = {
    '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4',
    '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
    'ₐ': 'a', 'ₑ': 'e', 'ₕ': 'h', 'ᵢ': 'i', 'ⱼ': 'j', 'ₖ': 'k',
    'ₗ': 'l', 'ₘ': 'm', 'ₙ': 'n', 'ₒ': 'o', 'ₚ': 'p', 'ᵣ': 'r',
    'ₛ': 's', 'ₜ': 't', 'ᵤ': 'u', 'ᵥ': 'v', 'ₓ': 'x',
}
SUP_MAP = {
    '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4',
    '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9',
    '⁺': '+', '⁻': '-', '⁼': '=', '⁽': '(', '⁾': ')',
    'ᵃ': 'a', 'ᵇ': 'b', 'ᶜ': 'c', 'ᵈ': 'd', 'ᵉ': 'e', 'ᶠ': 'f',
    'ᵍ': 'g', 'ʰ': 'h', 'ⁱ': 'i', 'ʲ': 'j', 'ᵏ': 'k', 'ˡ': 'l',
    'ᵐ': 'm', 'ⁿ': 'n', 'ᵒ': 'o', 'ᵖ': 'p', 'ʳ': 'r', 'ˢ': 's',
    'ᵗ': 't', 'ᵘ': 'u', 'ᵛ': 'v', 'ʷ': 'w', 'ˣ': 'x', 'ʸ': 'y', 'ᶻ': 'z',
    'ᴿ': '\\mathrm{R}', 'ᵀ': 'T',
}

CHAR_REPL = {
    'σ': r'\sigma ', 'Σ': r'\Sigma ',
    'ρ': r'\rho ', 'λ': r'\lambda ', 'δ': r'\delta ', 'τ': r'\tau ',
    'β': r'\beta ', 'φ': r'\phi ', 'ϕ': r'\varphi ',
    'ε': r'\epsilon ', 'κ': r'\kappa ', 'ω': r'\omega ',
    'θ': r'\theta ', 'π': r'\pi ', 'γ': r'\gamma ', 'ɣ': r'\gamma ',
    'Φ': r'\Phi ', 'Ψ': r'\Psi ', 'Ω': r'\Omega ', 'Λ': r'\Lambda ', 'Π': r'\Pi ', 'Θ': r'\Theta ', 'ψ': r'\psi ',
    'ℏ': r'\hbar ', '∂': r'\partial ', '∇': r'\gradv ',
    '∈': r'\in ', '∉': r'\notin ',
    '≤': r'\leq ', '≥': r'\geq ', '≠': r'\neq ', '≈': r'\approx ',
    '±': r'\pm ', '∓': r'\mp ',
    '×': r'\times ', '⨯': r'\times ', '•': r'\cdot ',
    '∑': r'\sum ', '∫': r'\int ', '∞': r'\infty ',
    '𝐈': r'\mathrm{I}',
    '𝓔': r'\mathbf E', '𝓑': r'\mathbf B', '𝓕': r'\mathcal{F}',
    'ℒ': r'\mathcal{L}', 'ℋ': r'\mathcal{H}',
    '□': r'\square ',
    '→': r'\Rightarrow ', '←': r'\leftarrow ', '⟨': r'\langle ', '⟩': r'\rangle ',
    'ℝ': r'\mathbb{R}', 'ℂ': r'\mathbb{C}',
    '½': r'\tfrac{1}{2}', '⅓': r'\tfrac{1}{3}', '¼': r'\tfrac{1}{4}', '⅛': r'\tfrac{1}{8}',
    '′': "'",
}

SIGMA_VECTOR_RE = re.compile(r'\\sigma(?!\s*_)')
SAFE_VECTOR_IDS = {'A', 'B', 'C', 'J', 'P', 'r', 'u', 'v'}
VECTOR_LHS_CANDIDATES = {'A', 'B', 'C', 'J', 'P', 'V', 'u', 'v', 'a', 'b', 'g', 'k', 'n', 'q', 'r'}
NONVECTOR_IDS = {
    'X', 'Y', 'Z', 'U',
    'i', 'scalar', 'vector', 'par', 'prp', 'real', 'img',
    'diag', 'det', 'exp', 'sin', 'cos', 'log', 'ln', 'tr',
}

# Keep this narrow so we really collect equation-bearing lines.
EQ_RE = re.compile(r"(=|:=|←|→|\bdet\(|\bdiag\(|\bexp\(|\bln\(|\blog\(|\bsin\(|\bcos\(|\bif\s+|\[.*\])")

PROSE_PREFIXES = (
    'using ',
    'with ',
    'for ',
    'if ',
    'note:',
    'let',
    'components',
    'real scalar:',
    'img scalar:',
    'real vector:',
    'img vector:',
    'dot product',
    'square vector notation',
    'scalar norm',
    'vector norm',
    'cross product',
    'connection to',
    'bivectors',
    'kronecker delta',
    'levi-civita symbol',
    'conjugate operation',
    'reverse operation',
    'projector products',
    'nilpotence',
    'action with projectors',
    'reconstruction',
    'eigenvalue decomposition',
    'matrix product',
    'spacetime',
    'complex conjugate',
    'complex conjugate transpose',
    'same-side projector',
    'cross-side projector',
    'same–side projector',
    'cross–side projector',
    'raising and lowering about z',
    '(ladder operators)',
    'rotation',
    'vectors',
)

COMPONENT_ROW_RE = re.compile(
    r'^(real scalar|img scalar|real vector|img vector|scalar|vector|real|img|bivector|bi-vector|pseudo-scalar|pseudoscalar)\s*:\s*(.*)$',
    re.IGNORECASE,
)

COMPONENT_BLOCK_COUNTER = 0
GRADE_TABLE_COUNTER = 0
COMPONENT_DISPLAY_LABELS = {
    'real scalar': 'real scalar',
    'img scalar': 'imaginary scalar',
    'real vector': 'real vector',
    'img vector': 'imaginary vector',
    'scalar': 'scalar',
    'vector': 'vector',
    'real': 'real part',
    'img': 'imaginary part',
    'bivector': 'bivector',
    'pseudo-scalar': 'pseudo-scalar',
}
GRADE_LABELS = {'scalar', 'vector', 'bivector', 'pseudo-scalar'}

QUALIFICATION_RULES = (
    (
        re.compile(r'^Robertson inequality$', re.IGNORECASE),
        r'This is not the standard Robertson or Robertson--Schr\"odinger uncertainty relation; the standard result is a state-variance inequality rather than the displayed paravector comparison.',
    ),
    (
        re.compile(r"^⟪\(du • σ\) \(du' • σ\) ⟫ ≥"),
        r'This displayed inequality should not be presented as the standard Robertson relation.',
    ),
    (
        re.compile(r"^ds² = dx' • dx'$"),
        r'This diagonalized form suppresses the Lorentzian signature; one must keep the indefinite sign structure explicit.',
    ),
    (
        re.compile(r'^T = ½ F Fᴿ$'),
        r'This is a frame-dependent energy-density/Poynting paravector, not by itself the full electromagnetic stress-energy tensor.',
    ),
    (
        re.compile(r'^\( W\(U\) – 𝐈 \) Ψ = 0'),
        r'This is an on-shell projector condition after normalization by mass, not the general Dirac equation.',
    ),
    (
        re.compile(r'^= a² –2 i \(b • σ\) \+ \(–i\)² \(b • σ\)²$'),
        r'The linear term is missing a factor of \(a\); it should be \(-2 i a (b\cdot\sigv)\).',
    ),
)

INTEGRATED_GROUPS = [
    {
        'title': 'Expanded Algebraic Foundations',
        'label': 'sec:integrated-algebra',
        'slug': 'algebra',
        'intro': (
            r'Section \ref{sec:algebra} introduced the compact algebraic toolkit used throughout the manuscript. '
            r'The present section unfolds that toolkit in more detail by recording the explicit commutator rules, '
            r'vector norms, involution tables, inverse formulas, and functional identities that stand behind '
            r'Eqs.\ \eqref{eq:basis}--\eqref{eq:functional-calculus}.'
        ),
        'blocks': [
            'commutator',
            'let',
            'σ algebra',
            'sigma conjugate operation',
            'reverse operation',
            'sign table',
            'σ vector form',
            'real products',
            'real exponent',
            'real logarithm',
            'paravectors: complex scalar + complex vector',
            'complex products',
            'paravector metric (not a norm)',
            'inverse',
            'minimal polynomial',
            'complex exponent',
            'complex fields',
            'functions of u • σ,  u² = 1',
            'functions of X',
        ],
    },
    {
        'title': 'Expanded Space-Time and Particle Kinematics',
        'label': 'sec:integrated-kinematics',
        'slug': 'kinematics',
        'intro': (
            r'Section \ref{sec:spacetime} established the paravector interpretation of space-time, proper velocity, '
            r'proper momentum, and free-particle dynamics. The following subsections write out the interval classes, '
            r'component decompositions, and variational formulas that sit behind Eqs.\ \eqref{eq:interval}--'
            r'\eqref{eq:free-lagrangian-low}.'
        ),
        'blocks': [
            'space–time numbers',
            'if ⟪dX⟫² > 0, the interval is "time-like"',
            '4 velocity  ← a.k.a  proper velocity',
            '4 acceleration ← a.k.a  proper acceleration',
            '4 momentum ← a.k.a  proper momentum',
            'Lagrangian',
            'action of a free particle =',
            'non-relativistic approximation',
        ],
    },
    {
        'title': 'Expanded Transformation Theory',
        'label': 'sec:integrated-transforms',
        'slug': 'transforms',
        'intro': (
            r'Section \ref{sec:transforms} presented the concise rotor and boost formulas. The following notes make the '
            r'transformation mechanism more explicit by expanding basis changes, Lorentz boosts, field transforms, and '
            r'composition laws associated with Eqs.\ \eqref{eq:general-transform}--\eqref{eq:field-transform}.'
        ),
        'blocks': [
            'transforms',
            'transformed basis',
            'Lorentz boost ( T = L )',
            'Lorentz transform of fields',
            'rotation ( T = R )',
            'compositions',
        ],
    },
    {
        'title': 'Expanded Electromagnetic Development',
        'label': 'sec:integrated-em',
        'slug': 'em',
        'intro': (
            r'Section \ref{sec:em} gave the compact STA formulation of electromagnetism. The following material expands '
            r'that formulation into explicit component identities for the field, wave operator, potentials, sources, '
            r'gauge transformations, Lorentz force, and action principles behind Eqs.\ \eqref{eq:dalembert}--'
            r'\eqref{eq:energy-momentum-paravector}.'
        ),
        'blocks': [
            'electricity & magnetism',
            'wave operator',
            'potential fields',
            'density from potential',
            'gauge transform',
            '4 force ← a.k.a  proper force',
            'Lagrangian derivation',
            'stress–energy of EM?',
        ],
    },
    {
        'title': 'Expanded Relativistic Quantum and Spin Development',
        'label': 'sec:integrated-qm',
        'slug': 'qm',
        'intro': (
            r'Sections \ref{sec:qm} and \ref{sec:spinrep} introduced the relativistic wave equations and their spinor '
            r'splittings in compact form. The following development writes out the detailed scalar and spinor '
            r'factorizations, root-operator interpretations, gauge couplings, low-energy limits, and spin projectors '
            r'used behind Eqs.\ \eqref{eq:kg-free}--\eqref{eq:energy-spin-projectors}.'
        ),
        'blocks': [
            'relativistic quantum mechanics',
            'spin 0 wave mechanics',
            'positive and negative frequency',
            'Interpreting root operators with Fourier',
            'spin ½ wave mechanics',
            'Klein–Gordon with potential fields',
            'Dirac with potential field',
            'gauge transform',
            'Dirac equation components',
            'low energy spin 0',
            'low energy spin ½',
            'spin - rotors',
            'spin - projectors',
            'spin uncertainty',
            'raising and lowering about z',
            'spin up and spin down',
            '4-state solution to Dirac equation',
        ],
    },
    {
        'title': 'Geometric Extensions and Representations',
        'label': 'sec:integrated-representations',
        'slug': 'representations',
        'intro': (
            r'The manuscript closes the integrated development with extensions that connect the STA presentation to '
            r'quaternions, generalized space-time forms, metric diagonalization attempts, and explicit matrix '
            r'representations. These notes are best read after the main physical development has already been absorbed.'
        ),
        'blocks': [
            'quaternions',
            'generalized spacetime',
            'metric',
            'matrix form',
        ],
    },
]

SKIP_INTEGRATED_TITLES = {
    'Initial Definitions',
    'SUMMARY OF RESULTS BELOW',
}


def escape_tex_text(s: str) -> str:
    s = s.replace('⟪', '<<').replace('⟫', '>>')
    s = s.replace('\\', r'\textbackslash{}')
    s = s.replace('&', r'\&').replace('%', r'\%').replace('#', r'\#')
    s = s.replace('_', r'\_').replace('{', r'\{').replace('}', r'\}')
    return s


def convert_subsup(s: str) -> str:
    out = []
    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch in SUB_MAP:
            j = i
            buf = []
            while j < n and s[j] in SUB_MAP:
                buf.append(SUB_MAP[s[j]])
                j += 1
            out.append('_{' + ''.join(buf) + '}')
            i = j
            continue
        if ch == '*' or ch in SUP_MAP:
            j = i
            buf = []
            while j < n and (s[j] == '*' or s[j] in SUP_MAP):
                if s[j] == '*':
                    buf.append('*')
                else:
                    buf.append(SUP_MAP[s[j]])
                j += 1
            out.append('^{' + ''.join(buf) + '}')
            i = j
            continue
        out.append(ch)
        i += 1
    return ''.join(out)


def replace_sqrt(s: str) -> str:
    out: list[str] = []
    i = 0
    n = len(s)
    while i < n:
        if s[i] != '√':
            out.append(s[i])
            i += 1
            continue

        if i + 1 < n and s[i + 1] == '(':
            depth = 1
            j = i + 2
            while j < n and depth > 0:
                if s[j] == '(':
                    depth += 1
                elif s[j] == ')':
                    depth -= 1
                j += 1
            if depth == 0:
                inner = replace_sqrt(s[i + 2:j - 1])
                out.append(r'\sqrt{' + inner + '}')
                i = j
                continue

        out.append(r'\sqrt{}')
        i += 1
    return ''.join(out)


def rhs_looks_math(s: str) -> bool:
    return bool(
        re.search(r'(=|:=|\d|σ|∂|Ψ|Φ|Ω|ρ|γ|ɣ|𝓔|𝓑|⟪|‖|\*|ᴿ|½|√)', s)
    )


def inline_rhs_looks_math(s: str) -> bool:
    return bool(
        re.search(r'(=|:=|\d|σ|∂|Ψ|Φ|Ω|ρ|γ|ɣ|𝓔|𝓑|⟪|‖|\*|ᴿ|½|√|ψ|•|×|₊|₋|ₛ|ᵥ|∈|ℝ|ℂ)', s)
    )


def replace_visible_sets(s: str) -> str:
    s = re.sub(r'\bdiag\s*\{([^{}]*)\}', lambda m: f"diag({m.group(1)})", s)
    return re.sub(
        r'([:=])\s*\{([^{}]*)\}',
        lambda m: f"{m.group(1)} \\left\\{{{m.group(2)}\\right\\}}",
        s,
    )


def vector_token_base(token: str) -> str:
    token = token.strip().replace(' ', '')
    token = re.sub(r"\([^)]*\)$", '', token)
    token = token.rstrip("'")
    if token.startswith('d') and len(token) > 1 and token[1:].isalpha():
        token = token[1:]
    return token


def collect_vector_ids(s: str) -> set[str]:
    vector_ids = set(SAFE_VECTOR_IDS)
    unbolded_prefix = r'(?<!\\mathbf )(?<!\\mathbf\{)(?<!\\bm\{)(?<!\\)'

    for pattern in (
        re.compile(unbolded_prefix + r"([A-Za-z]+(?:\([^)]*\))?(?:')?)\s*\\cdot\s*\\sigv"),
        re.compile(unbolded_prefix + r'([A-Za-z]+(?:\([^)]*\))?)\s*\\cdot\s*\\sigv'),
        re.compile(r"\\sigv\s*\\cdot\s*" + unbolded_prefix + r"([A-Za-z]+(?:\([^)]*\))?(?:')?)"),
        re.compile(r'\\sigv\s*\\cdot\s*' + unbolded_prefix + r'([A-Za-z]+(?:\([^)]*\))?)'),
        re.compile(unbolded_prefix + r"([A-Za-z]+(?:')?)\s*\\times\s*" + unbolded_prefix + r"([A-Za-z]+(?:')?)"),
        re.compile(unbolded_prefix + r'([A-Za-z]+)\s*\\times\s*' + unbolded_prefix + r'([A-Za-z]+)'),
        re.compile(unbolded_prefix + r"([A-Za-z]+(?:')?)\s*\\cdot\s*" + unbolded_prefix + r"([A-Za-z]+(?:')?)"),
        re.compile(unbolded_prefix + r'([A-Za-z]+)\s*\\cdot\s*' + unbolded_prefix + r'([A-Za-z]+)'),
        re.compile(r"\\lVert\s*" + unbolded_prefix + r"([A-Za-z]+(?:')?)\s*\\rVert"),
        re.compile(r'\\lVert\s*' + unbolded_prefix + r'([A-Za-z]+)\s*\\rVert'),
        re.compile(unbolded_prefix + r'([A-Za-z]+)\s*\\in\s*\\mathbb\{(?:R|C)\}\^\{?3\}?'),
        re.compile(r'\\partial\s*_\{?\s*[trv]\s*\}?\s*' + unbolded_prefix + r'([ABCPJabgkqrvu])'),
    ):
        for match in pattern.finditer(s):
            for group in match.groups():
                if group:
                    vector_ids.add(vector_token_base(group))

    lhs_match = re.match(r'^([A-Za-z]+)\s*[:=]+', s)
    if lhs_match:
        lhs = vector_token_base(lhs_match.group(1))
        rhs = s[lhs_match.end():]
        if lhs in VECTOR_LHS_CANDIDATES:
            if (
                re.search(r'\\mathbb\{(?:R|C)\}\^\{?3\}?', rhs)
                or r'\operatorname{vector}' in rhs
                or r'\lVert' in rhs
                or r'\times' in rhs
                or r'\sigv' in rhs
                or r'\partial _{r}' in rhs
                or r'\partial_{r}' in rhs
                or r'\partial _{v}' in rhs
                or r'\partial_{v}' in rhs
                or r'\delr' in rhs
                or r'\delv' in rhs
                or r'\nabla' in rhs
                or r'\gradv' in rhs
                or re.match(
                    rf'\s*\(\s*{re.escape(lhs)}_\{{?1\}}?\s*,\s*{re.escape(lhs)}_\{{?2\}}?\s*,\s*{re.escape(lhs)}_\{{?3\}}?\s*\)',
                    rhs,
                )
            ):
                vector_ids.add(lhs)

    if re.match(r'^q\s*[:=]+\s*-?i\s*\\sigv', s):
        vector_ids.add('q')

    return {ident for ident in vector_ids if ident and ident not in NONVECTOR_IDS}


def bold_identifier(s: str, ident: str) -> str:
    pattern = re.compile(
        rf'(?<!\\mathbf )(?<!\\mathbf\{{)(?<!\\bm\{{)(?<![A-Za-z]){re.escape(ident)}(?=(?:\'|\\prime)?(?![A-Za-z0-9_]))'
    )
    return pattern.sub(rf'\\mathbf {ident}', s)


def bold_vector_differentials(s: str, vector_ids: set[str]) -> str:
    for ident in sorted(vector_ids, key=len, reverse=True):
        pattern = re.compile(rf'(?<![A-Za-z\\])d{re.escape(ident)}(?=(?:\'|\\prime)?(?![A-Za-z_]))')
        s = pattern.sub(rf'd\\mathbf {ident}', s)
    return s


def bold_vector_derivative_subscripts(s: str, vector_ids: set[str]) -> str:
    for ident in sorted(vector_ids, key=len, reverse=True):
        if len(ident) != 1:
            continue
        pattern = re.compile(rf'\\partial\s*_\{{?\s*{re.escape(ident)}\s*\}}?')
        s = pattern.sub(rf'\\partial _{{\\mathbf {ident}}}', s)
    return s


def normalize_vector_notation(s: str) -> str:
    vector_ids = collect_vector_ids(s)
    s = bold_vector_differentials(s, vector_ids)
    s = bold_vector_derivative_subscripts(s, vector_ids)
    for ident in sorted(vector_ids, key=len, reverse=True):
        s = bold_identifier(s, ident)
    for old, new in (
        (r'\partial _{\mathbf r}', r'\delr{}'),
        (r'\partial_{\mathbf r}', r'\delr{}'),
        (r'\partial _{r}', r'\delr{}'),
        (r'\partial_{r}', r'\delr{}'),
        (r'\partial _{\mathbf v}', r'\delv{}'),
        (r'\partial_{\mathbf v}', r'\delv{}'),
        (r'\partial _{v}', r'\delv{}'),
        (r'\partial_{v}', r'\delv{}'),
        (r'\nabla', r'\gradv'),
    ):
        s = s.replace(old, new)
    s = re.sub(
        r'Q Q\^\{\\mathrm\{R\}\}\s*([&:]?=)\s*a\^\{2\}\s*\+\s*b\^\{2\}',
        lambda m: f"Q Q^{{\\mathrm{{R}}}} {m.group(1)} a^{{2}} + \\mathbf b^{{2}}",
        s,
    )
    return s


def romanize_equation_phrases(s: str) -> str:
    if (
        'symbolic vector' in s
        and '\\sigv' in s
        and '\\sigma _{1}' in s
        and '\\sigma _{2}' in s
        and '\\sigma _{3}' in s
    ):
        return r'\sigv = (\sigma _{1}, \sigma _{2}, \sigma _{3})\qquad \eqtext{vector of algebraic objects}.'

    s = re.sub(r'\\Rightarrow\s+commutative\b', r'\\Rightarrow \\text{commutative}', s)
    s = re.sub(r'\\qquad\s+different than\s+', r'\\qquad \\text{different from }', s)
    s = re.sub(r'\\text\{ if \}\s+any two indices are equal', r'\\text{ if any two indices are equal}', s)
    s = s.replace(r'non-commutative algebra \text{ if }', r'\text{non-commutative algebra if }')
    s = s.replace(r'constant \text{ if }', r'\text{constant if }')
    s = s.replace(r'\qquad 4 force', r'\qquad \text{4-force}')
    s = re.sub(r'(?<!\\text\{)3-element vector\b', r'\\text{3-element vector}', s)
    s = re.sub(r'(?<!\\text\{)real scalar \+ complex vector\b', r'\\text{real scalar + complex vector}', s)
    s = re.sub(r'(?<!\\text\{)real scalar \+ real vector\b', r'\\text{real scalar + real vector}', s)
    s = re.sub(r'(?<!\\text\{)total energy\b', r'\\text{total energy}', s)
    s = re.sub(r'(?<!\\text\{)kinetic \+ rest mass energy\b', r'\\text{kinetic + rest mass energy}', s)
    s = re.sub(r'(?m)^1\s*&=\s*vector,', r'1 &= \\text{vector},', s)
    s = re.sub(r'(?m)^2\s*&=\s*bi-vector,', r'2 &= \\text{bi-vector},', s)
    s = re.sub(r'(?m)^3\s*&=\s*pseudo-scalar,', r'3 &= \\text{pseudo-scalar},', s)
    s = re.sub(r'\bin\s+(?=\\mathbb)', lambda _m: r'\in ', s)
    s = s.replace(r'\text{complex scalar}', r'\eqtext{complex scalar}')
    s = s.replace(r'\text{real scalar}', r'\eqtext{real scalar}')
    s = s.replace(r'\text{complex vector}', r'\eqtext{complex vector}')
    s = s.replace(r'\text{real vector}', r'\eqtext{real vector}')
    s = s.replace(r'\text{vector of algebraic objects}', r'\eqtext{vector of algebraic objects}')
    s = s.replace(r'\text{3-element vector}', r'\eqtext{3-element vector}')
    s = s.replace(r'\text{real scalar + complex vector}', r'\eqtext{real scalar + complex vector}')
    s = s.replace(r'\text{real scalar + real vector}', r'\eqtext{real scalar + real vector}')
    s = s.replace(r'scalar^{2}', r'\eqtext{scalar}^{2}')
    s = s.replace(r'vector^{2}', r'\vectortext^{2}')
    s = s.replace(r'scalar + vector \cdot \sigv', r'\eqtext{scalar} + \vectortext \cdot \sigv')
    s = s.replace(r'scalar - vector \cdot \sigv', r'\eqtext{scalar} - \vectortext \cdot \sigv')
    return s


def convert_expr(line: str) -> str:
    s = line.strip()
    s = s.replace('–', '-').replace('—', '-').replace('−', '-')
    s = s.replace('“', '"').replace('”', '"').replace('’', "'").replace('‘', "'")
    if '←' in s:
        lhs, rhs = s.split('←', 1)
        lhs_tex = convert_expr(lhs.strip())
        rhs = rhs.strip()
        if not rhs:
            return lhs_tex
        if rhs_looks_math(rhs):
            return lhs_tex + r'\qquad ' + convert_expr(rhs)
        return lhs_tex + r'\qquad \text{' + escape_tex_text(rhs) + '}'

    s = replace_visible_sets(s)
    s = convert_subsup(s)
    s = replace_sqrt(s)

    # Norm and double-angle bracket handling.
    norm_open = True
    tmp = []
    for ch in s:
        if ch == '‖':
            tmp.append(r'\lVert ' if norm_open else r' \rVert')
            norm_open = not norm_open
        elif ch == '⟪':
            tmp.append(r'\langle\!\langle ')
        elif ch == '⟫':
            tmp.append(r' \rangle\!\rangle')
        else:
            tmp.append(ch)
    s = ''.join(tmp)

    for src, dst in CHAR_REPL.items():
        s = s.replace(src, dst)

    # Standalone sigma denotes the symbolic sigma-vector and should follow the
    # document's bold vector notation. Indexed basis elements \sigma_n remain
    # unbolded.
    s = SIGMA_VECTOR_RE.sub(r'\\sigv', s)
    s = normalize_vector_notation(s)

    bare_op_map = {
        'exp': r'\exp',
        'ln': r'\ln',
        'log': r'\log',
        'sin': r'\sin',
        'cos': r'\cos',
        'sinh': r'\sinh',
        'cosh': r'\cosh',
        'tanh': r'\tanh',
        'sech': r'\operatorname{sech}',
        'csch': r'\operatorname{csch}',
        'coth': r'\coth',
        'arcsin': r'\arcsin',
        'arccos': r'\arccos',
        'arctan': r'\arctan',
        'arcsinh': r'\operatorname{arcsinh}',
        'arccosh': r'\operatorname{arccosh}',
        'arctanh': r'\operatorname{arctanh}',
    }
    for src, dst in sorted(bare_op_map.items(), key=lambda item: -len(item[0])):
        s = re.sub(rf'(?<![A-Za-z\\]){src}(?![A-Za-z])', lambda _m, repl=dst: repl, s)

    call_op_map = {
        'diag': r'\operatorname{diag}',
        'det': r'\operatorname{det}',
        'real': r'\eqtext{real}',
        'img': r'\eqtext{img}',
        'scalar': r'\eqtext{scalar}',
        'vector': r'\eqtext{vector}',
        'par': r'\eqtext{par}',
        'prp': r'\eqtext{prp}',
        'tr': r'\eqtext{tr}',
    }
    for src, dst in call_op_map.items():
        s = re.sub(rf'\b{src}\s*\(', lambda _m, repl=dst: repl + '(', s)

    s = re.sub(r'(?<=\s)if and only if(?=\s)', r'\\text{ if and only if }', s)
    s = re.sub(r'(?<=\s)if(?=\s)', r'\\text{ if }', s)
    s = re.sub(r'(?<=\s)else(?=\s)', r'\\text{ else }', s)
    s = s.replace('to commutative', r'\Rightarrow \text{commutative}')

    # Light cleanup of repeated spaces.
    s = re.sub(r'\s+', ' ', s).strip()

    s = romanize_equation_phrases(s)

    # Escape TeX-sensitive chars that may appear in raw lines.
    s = s.replace('&', r'\&').replace('%', r'\%').replace('#', r'\#')

    return s


def is_equation_line(line: str) -> bool:
    s = line.strip()
    if not s:
        return False
    # Strong signal: equation or relation symbols.
    if EQ_RE.search(s):
        return True
    # Additional symbol-based fallback for symbolic lines.
    for token in ('σ', '∂', 'Ψ', 'Φ', 'Ω', 'ρ', '𝓔', '𝓑', '⟪', '‖', 'ɣ'):
        if token in s:
            return True
    return False


def is_display_math_line(line: str) -> bool:
    s = line.strip()
    if not is_equation_line(s):
        return False
    lower = s.lower()
    if any(lower.startswith(prefix) for prefix in PROSE_PREFIXES):
        return False
    if s.endswith(':'):
        return False
    return True


def qualification_note(line: str) -> str | None:
    for pattern, note in QUALIFICATION_RULES:
        if pattern.search(line.strip()):
            return note
    return None


def integrated_display_title(title: str, group_title: str) -> str:
    title_map = {
        'commutator': 'Commutator Algebra',
        'let': 'Vector Notation, Norms, and Cross Products',
        'σ algebra': 'Sigma Algebra',
        'sigma conjugate operation': 'Conjugation of Sigma Products',
        'reverse operation': 'Reversion of Sigma Products',
        'sign table': 'Sign Table for the Involutions',
        'σ vector form': 'Vector Form and Explicit Sigma Products',
        'real products': 'Products of Real Scalars, Vectors, and Paravectors',
        'real exponent': 'Exponential of a Real Vector',
        'real logarithm': 'Logarithm of a Real Paravector',
        'paravectors: complex scalar + complex vector': 'Complex Paravectors',
        'complex products': 'Products of Complex Paravectors',
        'paravector metric (not a norm)': 'Paravector Metric and Interval Form',
        'inverse': 'Inverse Formulas',
        'minimal polynomial': 'Minimal Polynomial and Spectral Decomposition',
        'complex exponent': 'Exponentials of Complex Paravectors',
        'complex fields': 'Complex Field Decomposition',
        'space–time numbers': 'Spacetime Paravectors',
        'if ⟪dX⟫² > 0, the interval is "time-like"': 'Timelike, Null, and Spacelike Classification',
        '4 velocity  ← a.k.a  proper velocity': 'Proper Velocity',
        '4 acceleration ← a.k.a  proper acceleration': 'Proper Acceleration',
        '4 momentum ← a.k.a  proper momentum': 'Proper Momentum',
        'Lagrangian': 'Relativistic Lagrangian',
        'action of a free particle =': 'Free-Particle Action',
        'non-relativistic approximation': 'Nonrelativistic Expansion',
        'electricity & magnetism': 'Electric and Magnetic Fields',
        'wave operator': 'Wave Operator and d’Alembertian',
        'potential fields': 'Potential Fields',
        'density from potential': 'Sources Derived from the Potential',
        '4 force ← a.k.a  proper force': 'Proper Force and Lorentz Dynamics',
        'Lagrangian derivation': 'Variational Derivation of the Lorentz Force',
        'stress–energy of EM?': 'Electromagnetic Energy-Momentum',
        'relativistic quantum mechanics': 'Relativistic Quantum-Mechanical Setting',
        'spin 0 wave mechanics': 'Spin-0 Wave Mechanics',
        'positive and negative frequency': 'Positive and Negative Frequency Sectors',
        'Interpreting root operators with Fourier': 'Root Operators and Fourier Interpretation',
        'spin ½ wave mechanics': 'Spin-1/2 Wave Mechanics',
        'Klein–Gordon with potential fields': 'Klein-Gordon Equation with Potentials',
        'Dirac with potential field': 'Dirac Equation with Potentials',
        'Dirac equation components': 'Dirac Equation in Components',
        'low energy spin 0': 'Low-Energy Spin-0 Limit',
        'low energy spin ½': 'Low-Energy Spin-1/2 Limit',
        'spin - rotors': 'Spin via Rotors',
        'spin - projectors': 'Spin Projectors',
        'spin uncertainty': 'Spin Uncertainty and Commutator Bounds',
        'raising and lowering about z': 'Raising and Lowering About the z Axis',
        'spin up and spin down': 'Spin-Up and Spin-Down Decomposition',
        '4-state solution to Dirac equation': 'Four-Sector Decomposition of the Dirac System',
        'transforms': 'General Transformation Identities',
        'transformed basis': 'Transformed Basis Elements',
        'Lorentz boost ( T = L )': 'Lorentz Boosts',
        'Lorentz transform of fields': 'Lorentz Transformation of the Field',
        'rotation ( T = R )': 'Spatial Rotations',
        'compositions': 'Compositions of Boosts and Rotations',
        'quaternions': 'Quaternion Correspondence',
        'functions of u • σ,  u² = 1': 'Functions of a Unit-Direction Element',
        'functions of X': 'Functions of a Paravector',
        'generalized spacetime': 'Generalized Spacetime Forms',
        'metric': 'Metric and Diagonalization Notes',
        'matrix form': 'Matrix Realization',
    }
    if title == 'gauge transform':
        if 'Electromagnetic' in group_title:
            return 'Gauge Transform of the Potential'
        return 'Gauge Covariance in the Dirac System'
    return title_map.get(title, title)


def integrated_ref_text(title: str, group_title: str) -> str:
    ref_map = {
        'commutator': r'Eqs.\ \eqref{eq:sigma-product} and \eqref{eq:dot-cross}',
        'σ algebra': r'Eqs.\ \eqref{eq:basis} and \eqref{eq:sigma-product}',
        'sigma conjugate operation': r'Eqs.\ \eqref{eq:sigma-involutions} and \eqref{eq:i-involutions}',
        'reverse operation': r'Eqs.\ \eqref{eq:sigma-involutions} and \eqref{eq:i-involutions}',
        'sign table': r'Eqs.\ \eqref{eq:sigma-involutions} and \eqref{eq:i-involutions}',
        'σ vector form': r'Eq.\ \eqref{eq:dot-cross}',
        'paravectors: complex scalar + complex vector': r'Eqs.\ \eqref{eq:paravector} and \eqref{eq:paravector-parts}',
        'complex products': r'Eq.\ \eqref{eq:paravector-product}',
        'paravector metric (not a norm)': r'Eqs.\ \eqref{eq:euclidean-form} and \eqref{eq:minkowski-form}',
        'inverse': r'Eq.\ \eqref{eq:minkowski-form}',
        'minimal polynomial': r'Eq.\ \eqref{eq:functional-calculus}',
        'complex exponent': r'Eqs.\ \eqref{eq:exp-vector} and \eqref{eq:exp-bivector}',
        'space–time numbers': r'Eqs.\ \eqref{eq:dX} and \eqref{eq:interval}',
        '4 velocity  ← a.k.a  proper velocity': r'Eq.\ \eqref{eq:proper-velocity}',
        '4 acceleration ← a.k.a  proper acceleration': r'Eq.\ \eqref{eq:proper-acceleration}',
        '4 momentum ← a.k.a  proper momentum': r'Eq.\ \eqref{eq:proper-momentum}',
        'Lagrangian': r'Eqs.\ \eqref{eq:free-lagrangian} and \eqref{eq:free-lagrangian-low}',
        'electricity & magnetism': r'Eqs.\ \eqref{eq:dalembert}, \eqref{eq:F-Phi}, and \eqref{eq:maxwell}',
        'wave operator': r'Eq.\ \eqref{eq:dalembert}',
        'potential fields': r'Eqs.\ \eqref{eq:potential-gradient} and \eqref{eq:E-B-from-potential}',
        'density from potential': r'Eq.\ \eqref{eq:maxwell}',
        '4 force ← a.k.a  proper force': r'Eq.\ \eqref{eq:lorentz-force}',
        'Lagrangian derivation': r'Eq.\ \eqref{eq:em-lagrangian}',
        'stress–energy of EM?': r'Eq.\ \eqref{eq:energy-momentum-paravector}',
        'spin 0 wave mechanics': r'Eq.\ \eqref{eq:kg-free}',
        'positive and negative frequency': r'Eq.\ \eqref{eq:frequency-split}',
        'spin ½ wave mechanics': r'Eqs.\ \eqref{eq:W-map} and \eqref{eq:dirac-free}',
        'Klein–Gordon with potential fields': r'Eqs.\ \eqref{eq:kg-coupled} and \eqref{eq:kg-expanded}',
        'Dirac with potential field': r'Eq.\ \eqref{eq:dirac-coupled}',
        'Dirac equation components': r'Eq.\ \eqref{eq:dirac-square}',
        'low energy spin 0': r'Eq.\ \eqref{eq:kg-coupled}',
        'low energy spin ½': r'Eqs.\ \eqref{eq:pauli-operator} and \eqref{eq:pauli}',
        'spin - rotors': r'Eq.\ \eqref{eq:rotor}',
        'spin - projectors': r'Eq.\ \eqref{eq:projectors}',
        'spin up and spin down': r'Eqs.\ \eqref{eq:pauli-spin-decomposition} and \eqref{eq:pauli-spin-split}',
        '4-state solution to Dirac equation': r'Eqs.\ \eqref{eq:energy-spin-projectors} and \eqref{eq:normalized-dirac}',
        'transforms': r'Eq.\ \eqref{eq:general-transform}',
        'Lorentz boost ( T = L )': r'Eq.\ \eqref{eq:boost}',
        'Lorentz transform of fields': r'Eq.\ \eqref{eq:field-transform}',
        'rotation ( T = R )': r'Eq.\ \eqref{eq:rotor}',
        'quaternions': r'Eqs.\ \eqref{eq:quaternion-basis} and \eqref{eq:quaternion-square}',
        'functions of u • σ,  u² = 1': r'Eqs.\ \eqref{eq:projectors} and \eqref{eq:functional-calculus}',
        'matrix form': r'Eqs.\ \eqref{eq:pauli-matrices} and \eqref{eq:matrix-X}',
    }
    if title == 'gauge transform':
        if 'Electromagnetic' in group_title:
            return r'Eq.\ \eqref{eq:gauge}'
        return r'Eq.\ \eqref{eq:dirac-gauge}'
    return ref_map.get(title, '')


def integrated_block_intro(title: str, group_title: str) -> str:
    display = integrated_display_title(title, group_title)
    ref_text = integrated_ref_text(title, group_title)
    custom = {
        'commutator': 'Before specializing to particular STA objects, it is helpful to record the algebraic behavior of the commutator itself, since many later identities are compact ways of saying that certain products commute or fail to commute in controlled ways.',
        'let': 'The elementary vector conventions are made fully explicit before introducing the sigma basis. That choice is retained here so that later products and norms can be read against familiar Euclidean notation.',
        'σ algebra': 'The sigma generators and their multiplication rules are written out here in full so that every later paravector identity can be traced back to explicit algebraic relations.',
        'products': 'This subsection expands the basic scalar, vector, and paravector products so that later identities can be read as consequences of one common multiplication rule rather than as isolated formulas.',
        'space–time numbers': 'At this point the paravector formalism is reinterpreted as spacetime. The following identities make that reinterpretation explicit by separating temporal and spatial parts and identifying the interval carried by the same algebraic product.',
        'if ⟪dX⟫² > 0, the interval is "time-like"': 'Once the interval has been written as a Minkowski-type quadratic form, its sign can be used to classify displacements. The following notes record that classification in the notation of the manuscript.',
        'electricity & magnetism': 'The next step is to spell out the electromagnetic field content that was previously presented in compact STA form. These identities make the relation between the vector fields and the single field paravector fully explicit.',
        'wave operator': "The wave operator is written out here in full so that its factorization into the d'Alembertian and the component field equations can be read directly.",
        'potential fields': 'The potential-based formulation is expanded here because it is the point at which gauge freedom, source reconstruction, and the emergence of the electric and magnetic fields all become transparent at once.',
        'minimal polynomial': 'The spectral data of a complex paravector are collected here because the later functional calculus, projector formulas, and exponential identities all depend on this algebraic factorization.',
        'positive and negative frequency': 'The frequency decomposition is one of the central interpretive steps in the scalar wave theory. The following formulas write out the split in a way that makes later Fourier and root-operator remarks easier to follow.',
        'Interpreting root operators with Fourier': 'Square-root operators can look formal if they are introduced without analytic context. This subsection explains them by switching to Fourier space, where their action reduces to multiplication by an explicit dispersion function.',
        'Klein–Gordon with potential fields': 'This subsection expands the minimally coupled Klein--Gordon equation in the presence of electromagnetic potentials. The point is to show explicitly how the compact operator form reproduces the familiar scalar equation with gauge-coupled derivatives and field-dependent commutators.',
        'Dirac equation components': 'The block Dirac equation becomes most transparent when its component equations are written out explicitly. Doing so makes it easier to compare the STA formulation with the familiar coupled two-component form.',
        'spin uncertainty': 'The commutator identities relevant to the uncertainty discussion are collected here in the same notation used elsewhere in the manuscript. The displayed inequality itself must still be read together with the red caution note that follows.',
        'metric': 'For comparison with standard linear-algebra notation, the paravector quadratic form is also written below in matrix language. The final diagonalization claim must still be read together with the red qualification that follows.',
        'quaternions': 'The quaternion correspondence is written out explicitly here to show how it sits inside the same multiplication table already developed for the sigma algebra.',
    }
    intro = custom.get(title)
    if intro is None:
        intro = f'The formulas below make {display.lower()} explicit and collect identities used later in the manuscript.'
    if ref_text:
        intro += f' See {ref_text}.'
    local_type_notes = {
        'transforms': r'Here \(T\) denotes a generic invertible transform element. The specific cases used later are \(R\) for spatial rotations and \(L\) for Lorentz boosts.',
        'commutator': r'Here \(X,Y\in\C\oplus\C^3\) are generic complex paravectors, and \([X,Y]:=XY-YX\) is their commutator.',
        'paravector metric (not a norm)': r'Here \(X,Y\in\R\oplus\R^3\) and \(Z\in\C\oplus\C^3\), so the displayed quadratic forms can be read against the real and complex paravector decompositions already established in the main text.',
        'space–time numbers': r"Here \(dt,dt',ds,s,\gamma,m\in\R\), \(d\mathbf r,d\mathbf r',\mathbf r,\mathbf v,\mathbf a,\mathbf P,\mathbf f\in\R^3\), and \(dX,U,A,mU\in\R\oplus\R^3\).",
        'if ⟪dX⟫² > 0, the interval is "time-like"': r"Here \(dX,dX'\in\R\oplus\R^3\), while \(dt,dt',ds,\tau\in\R\) and \(d\mathbf r,d\mathbf r'\in\R^3\).",
        '4 velocity  ← a.k.a  proper velocity': r'Here \(\mathbf v\in\R^3\), \(\gamma\in\R\), and \(U=dX/ds\in\R\oplus\R^3\) is the proper-velocity paravector.',
        '4 acceleration ← a.k.a  proper acceleration': r'Here \(\mathbf a\in\R^3\), \(s\in\R\), and \(A=dU/ds\in\R\oplus\R^3\) is the proper-acceleration paravector.',
        'Lorentz boost ( T = L )': r'Here \(\theta\in\R\), \(\mathbf u\in\R^3\) with \(\|\mathbf u\|=1\), and \(L\) is the boost element acting by conjugation on paravectors.',
        'rotation ( T = R )': r'Here \(\theta\in\R\), \(\mathbf u\in\R^3\) with \(\|\mathbf u\|=1\), and \(R\) is the rotor element acting by conjugation on paravectors.',
        'Lorentz transform of fields': r'Here \(F\in\{0\}\oplus\C^3\subset\C\oplus\C^3\), while \(\mathbf E,\mathbf B,\mathbf u,\mathbf v\in\R^3\) and \(\gamma,\theta\in\R\).',
        'electricity & magnetism': r'Here \(V,\rho,S,q\in\R\), \(\mathbf A,\mathbf E,\mathbf B,\mathbf J\in\R^3\), \(\Phi=V+\mathbf A\cdot\sigv\in\R\oplus\R^3\), and \(F=(\mathbf E+i\mathbf B)\cdot\sigv\in\{0\}\oplus\C^3\).',
        'potential fields': r'Here \(V,S,\rho,\lambda\in\R\), \(\mathbf A,\mathbf E,\mathbf B,\mathbf J\in\R^3\), \(\Phi\in\R\oplus\R^3\), and \(F\in\{0\}\oplus\C^3\subset\C\oplus\C^3\).',
        '4 force ← a.k.a  proper force': r'Here \(q,m\in\R\), \(F\in\{0\}\oplus\C^3\subset\C\oplus\C^3\), \(U\in\R\oplus\R^3\), and \(\mathbf E,\mathbf B,\mathbf v,\mathbf f\in\R^3\).',
        'Lagrangian derivation': r'Here \(X,U,P\in\R\oplus\R^3\), \(s,t,q,m\in\R\), \(\mathbf r,\mathbf v,\mathbf A\in\R^3\), and \(V\in\R\).',
        'relativistic quantum-mechanical setting': r'Here \(\psi=\psi_s+\psi_{\mathbf v}\cdot\sigv\in\C\oplus\C^3\) denotes a generic complex paravector state, while \(M\) is an operator acting on that state.',
        'spin 0 wave mechanics': r'Here \(\psi\in\C\) is a complex scalar field, \(m\in\R\), and \(\partial\) is the paravector differential operator.',
        'positive and negative frequency': r'Here \(\psi,\psi(\pm),c(\pm)\in\C\), while the Fourier parameters \(\omega_n\in\R\) and \(\mathbf k_n\in\R^3\) appear when the root operator is analyzed mode by mode.',
        'Interpreting root operators with Fourier': r'In the Fourier discussion, \(\psi\in\C\), the coefficients \(c_n\in\C\), the frequencies \(\omega_n\in\R\), and the wave vectors \(\mathbf k_n\in\R^3\).',
        'spin ½ wave mechanics': r'Here \(\bPsi=\{\bPsi_1;\bPsi_2\}\) is a two-component block field, \(Z\in\C\oplus\C^3\), and \(\Wmat(Z)\) is a \(2\times2\) block matrix with paravector entries.',
        'Klein–Gordon with potential fields': r'Here \(\psi\in\C\), \(\Phi=V+\mathbf A\cdot\sigv\in\R\oplus\R^3\), \(F\in\{0\}\oplus\C^3\subset\C\oplus\C^3\), \(V,S,q,m\in\R\), and \(\mathbf A,\mathbf E,\mathbf B\in\R^3\).',
        'Dirac with potential field': r'Here \(\bPsi=\{\bPsi_1;\bPsi_2\}\) is a two-component block field, \(\Phi\in\R\oplus\R^3\), and \(\Dmat\) is a block operator built from paravector entries.',
        'Dirac equation components': r'Here \(\psi=\psi_s+\psi_{\mathbf v}\cdot\sigv\in\C\oplus\C^3\), \(\Phi\in\R\oplus\R^3\), and \(F=F_{\mathbf v}\cdot\sigv\in\{0\}\oplus\C^3\) with \(F_{\mathbf v}\in\C^3\).',
        'low energy spin 0': r'Here \(\psi\in\C\), while \(E\) and \(\mathbf P\) remain the scalar and vector differential operators defined in the main text.',
        'low energy spin ½': r'Here \(\bPsi_1\) and \(\bPsi_2\) denote the block components of the Dirac field, and the low-energy reduction keeps the same real fields \(V\in\R\) and \(\mathbf A,\mathbf B\in\R^3\).',
        'spin - rotors': r'Here \(\mathbf u,\mathbf n\in\R^3\) are unit directions, and \(R\) is a rotor algebra element.',
        'spin - projectors': r'Here \(\mathbf u,\mathbf u^{\prime}\in\R^3\) are unit directions and \(\Pi(\pm)\) are the corresponding algebraic projectors.',
        'raising and lowering about z': r'Here \(\mathbf u=e_3\in\R^3\), and \(a(\pm)\) are algebra elements built from \(\sigma_1\) and \(\sigma_2\).',
        'spin up and spin down': r'Here \(\mathbf B,\mathbf u\in\R^3\), \(\|\mathbf u\|=1\), and \(\bPsi(\pm)\) are the projector-selected spin sectors of \(\bPsi\).',
        '4-state solution to Dirac equation': r'Here \(\epsilon,s\in\{+,-\}\), \(\Lmat(\epsilon)\) and \(\Pi(s)\) are the energy and spin projectors, and \(\bPsi(\epsilon,s)\) is the corresponding block-field sector.',
        'quaternions': r'Here \(a\in\R\), \(\mathbf b\in\R^3\), and \(Q=a+\mathbf b\cdot\mathbf q=a-i\mathbf b\cdot\sigv\in\C\oplus\C^3\) is the quaternion-paravector correspondence used in this section.',
        'metric': r'Here \(X\in\R\oplus\R^3\), \(dx\in\R\oplus\R^3\), and \(\Qmat(x)\) is the matrix-valued bilinear form used for comparison with the intrinsic paravector metric.',
    }
    type_note = local_type_notes.get(title)
    if type_note:
        intro += f' {type_note}'
    return intro


def consume_block(blocks: list[tuple[str, list[str]]], title: str) -> tuple[str, list[str]] | None:
    for idx, block in enumerate(blocks):
        if block[0] == title:
            return blocks.pop(idx)
    return None


def parse_blocks(lines: list[str]) -> list[tuple[str, list[tuple[int, str]]]]:
    blocks: list[tuple[str, list[tuple[int, str]]]] = []
    current_title = 'Initial Definitions'
    current_items: list[tuple[int, str]] = []
    expecting_title = False

    def flush() -> None:
        nonlocal current_items, current_title
        if current_items:
            blocks.append((current_title, current_items))
            current_items = []

    for idx, raw in enumerate(lines, start=1):
        line = raw.rstrip('\n')
        stripped = line.strip()

        if re.fullmatch(r'_+', stripped):
            flush()
            expecting_title = True
            continue

        if expecting_title and stripped:
            current_title = stripped
            expecting_title = False
            continue

        if is_equation_line(stripped):
            current_items.append((idx, stripped))

    flush()
    return blocks


def indent_tex(raw_line: str) -> str:
    leading = len(raw_line) - len(raw_line.lstrip(' '))
    if leading <= 0:
        return ''
    return rf'\hspace*{{{leading * 0.5:.1f}em}}'


def normalize_plain_text(s: str) -> str:
    s = s.replace('–', '-').replace('—', '-').replace('−', '-')
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def declaration_phrase(desc: str) -> str:
    mapped = {
        '3-element vector': 'a real three-vector',
        'scalar charge density': 'the scalar charge density',
        'vector current density': 'the current-density three-vector',
        'electric field': 'the electric field',
        'magnetic field': 'the magnetic field',
        'electro-magnetic field': 'the complex electromagnetic field',
        'algebraic object': 'a basis algebra element',
        'symbolic vector': 'the symbolic sigma-vector',
        'vector of bivectors': 'the vector of bivectors',
        'scalar operator': 'a scalar operator',
        'vector operator': 'a vector differential operator',
        'paravector': 'a paravector',
        'constant': 'a constant',
    }
    return mapped.get(desc.lower(), desc.lower())


def inline_math(expr: str) -> str:
    return r'\(' + convert_expr(expr) + r'\)'


def ensure_sentence(text: str) -> str:
    text = text.strip()
    if text.endswith(('.', '!', '?')):
        return text
    return text + '.'


def format_prose_line(stripped: str) -> str | None:
    s = normalize_plain_text(stripped)
    lower = s.lower()

    cue_map = {
        'let': 'Let the following notation be fixed.',
        'dot product': 'The dot product is defined as follows.',
        'addition': 'Addition is handled componentwise.',
        'multiplication': 'The corresponding multiplication law is',
        'cross product': 'The cross product is defined by',
        'square vector notation': 'We denote the squared vector magnitude by',
        'scalar norm': 'The scalar norm is defined by',
        'vector norm': 'The vector norm is defined by',
        'vectors': 'Introduce the following vectors.',
        'vector': 'Introduce the following vector notation.',
        'density paravector': 'Introduce the source-density paravector.',
        'this expands to': 'Expanding the left-hand side gives',
        'target relationships': 'The target component identities are',
        'consider': 'Consider the following identity.',
        'without s': 'Equivalently, without introducing \(S\), one may write',
        'define': 'Define the following object.',
        'summary': 'Collecting the definitions just obtained, one has',
        'components': 'Decomposing the expression into scalar and vector parts gives',
        'paravector metric': 'The quadratic paravector form is',
        'rotation': 'The corresponding rotation statement is',
        'wave equation': 'The corresponding wave equation is',
        "d'alembertian": "The d'Alembertian is",
        'd’alembertian': "The d'Alembertian is",
        'hamiltonian': 'The Hamiltonian is',
        'canonical momentum': 'The canonical momentum is',
        'relativistic energy': 'The relativistic energy is',
        'relativistic momentum': 'The relativistic momentum is',
        'relativistic velocity': 'The relativistic velocity is',
        'geometric velocity': 'The geometric velocity is',
        'geometric lorentz factor': 'The geometric Lorentz factor is',
        'geometric acceleration': 'The geometric acceleration is',
        'relativistic acceleration': 'The relativistic acceleration is',
        'bra-ket notation': 'Use the following bra-ket notation.',
        'expectation': 'The expectation value is',
        'normalized': 'Assume the state is normalized.',
        'mass conjugate operator': 'Introduce also the mass-conjugate operator.',
        'gives': 'The resulting relations are',
        'with': 'With these definitions fixed, one obtains',
        'projection operator': 'Define the projector family by',
        'projection operators': 'Define the projection operators by',
        'geometric velocity': 'The geometric velocity is',
        'unit momentum': 'Introduce the unit momentum paravector.',
        'unit momentum operator': 'Introduce the normalized momentum operator.',
        'on mass shell': 'On the mass shell one has',
        'solution': 'The resulting solution condition is',
        'add zero': 'Add a vanishing term to expose the desired decomposition.',
        'vector calc identity': 'Use the following vector-calculus identity.',
        'in one line': 'These identities can be compressed into the single relation',
        'in one formula': 'These terms combine into the single relation',
        'relativistic work': 'The time component gives the relativistic work law',
        'relativistic lorentz force law': 'The spatial component gives the relativistic Lorentz force law',
        'spacetime': 'The spacetime interpretation is',
        'complex conjugate': 'The complex-conjugation rules are',
        'complex conjugate transpose': 'The conjugate-transpose rules are',
        'reverse symmetry': 'The reverse symmetry can be written as',
        'differentials': 'Differentiating the transformed quantities gives',
        'eigenvalue decomposition': 'The eigenvalue decomposition is',
        'projection operators': 'Define the projection operators by',
        'projection operator': 'Define the projector by',
        'eigenvalues': 'The eigenvalues are',
        'eigen-projectors': 'The corresponding projectors are',
        'eigenvalue equation': 'The eigenvalue equation is',
        'projector decomposition': 'The projector decomposition is',
        'minimal polynomial': 'The minimal polynomial is',
        'properties': 'These operators satisfy the following properties.',
        'momentum': 'The momentum is',
        'zeeman effect': 'For a uniform static magnetic field, the Zeeman splitting becomes',
        'energy projector': 'Define the energy projectors by',
        'spin projector': 'Define the spin projectors by',
        're-arranging the last two': 'Rearranging the vector equations gives',
        'recovering elements': 'The original scalar, vector, real, and imaginary parts can now be recovered from these involutions as follows.',
        'where': 'From these identities, the resulting parameter formulas are',
        'dirac with potential field': 'Begin with the Dirac equation in the presence of a potential field.',
        'dirac operator with potential': 'Introduce the Dirac operator in the presence of a potential field.',
        '2 equations': 'This produces two coupled equations.',
        '2 equations without rest mass': 'After removing the rest phase, the system becomes',
        'peel off rest mass': 'Next factor out the rest-mass oscillation.',
        'plug it into second equation': 'Substituting this into the second equation yields',
        'plug into second equation': 'Substituting this into the second equation yields',
        'plug into the second equation': 'Substituting this into the second equation yields',
        'and multiply by m': 'After multiplying through by \(m\), one obtains',
        'exact solution': 'Keeping the full operator structure gives the exact equation',
        'low energy approx': 'At leading order in the low-energy expansion, this reduces to',
        'low energy approximation': 'At leading order in the low-energy expansion, this reduces to',
        'rest mass oscillation': 'Factor out the rest-mass oscillation.',
        'schrödinger equation': "The resulting Schr\"odinger equation is",
        'pauli equation': 'The resulting Pauli equation is',
        'klein–gordon equation': 'The corresponding Klein--Gordon equation is',
        'rodrigues’ formula': "Applying the rotor action gives Rodrigues' formula,",
        'unit vectors': 'Choose unit vectors as follows.',
        'with conjugate': 'Applying conjugation gives the companion relation',
        'recall': 'Recall the previously established definitions.',
        'kronecker delta': 'The Kronecker delta is defined by',
        'levi-civita symbol': 'The Levi-Civita symbol is defined by',
        'connection to imaginary numbers': 'The pseudoscalar connects the sigma algebra to the ordinary imaginary unit.',
        'bivectors': 'The bivector sector is given by',
        'all product components can be reduced to': 'Every product can be reduced to the following grade components.',
    }
    if lower in cue_map:
        return rf'\noindent\textit{{{ensure_sentence(cue_map[lower])}}}\par'

    if lower == '≠ □':
        return rf'\noindent\textit{{This differs from }}{inline_math("□")}\textit{{.}}\par'

    if lower == 'note:':
        return r'\noindent\emph{Note.} The following supporting identities are used in the surrounding derivation.\par'

    if lower.startswith('note:'):
        body = s.split(':', 1)[1].strip()
        if not body:
            return None
        return rf'\noindent\emph{{Note.}} {ensure_sentence(escape_tex_text(body))}\par'

    if lower in {'(measured along a path).', '(measured along a worldline).'}:
        return r'\noindent\emph{Here the quantity is measured along the worldline.}\par'

    if lower == 'e.g., non-inertial frame.':
        return r'\noindent\emph{This interpretation also applies in non-inertial settings.}\par'

    if lower == 'if boost or rotor,.':
        return r'\noindent\textit{For rotors and boosts, the inverse simplifies to}\par'

    if lower == 'conjugate and reverse.':
        return r'\noindent\textit{Applying conjugation and reversion gives}\par'

    if lower == 'energy & momentum.':
        return r'\noindent\textit{The transformed energy and momentum components are}\par'

    if lower == 'inertial motion:':
        return r'\noindent\textit{For inertial motion, one has}\par'

    if lower in {'some use the (-,+,+,+) signature', 'some use the  (-,+,+,+) signature'}:
        return r'\noindent\emph{Note.} Some authors instead use the metric signature \((- ,+,+,+)\).\par'

    if lower == 'some use □²':
        return r'\noindent\emph{Note.} Some authors denote the wave operator with a different box-symbol convention.\par'

    if lower in {'with σ-algebra.', 'with sigma-algebra.'}:
        return r'\noindent\emph{Note.} This formula is written using the sigma-algebra conventions fixed earlier in the manuscript.\par'

    if lower == 'to remain time-like, v² < 1.':
        return r'\noindent\emph{Note.} The timelike condition requires \(\mathbf v^{2}<1\).\par'

    if lower == 'the arc length s equals the proper time τ':
        return (
            r'\noindent The arc length '
            + inline_math("s")
            + r' agrees with the proper time '
            + inline_math(r'\tau')
            + r', so the same scalar parameter measures elapsed time along the worldline.\par'
        )

    if lower == 'dx⁻¹ does not exist' or lower == 'dx^-1 does not exist' or lower == 'dx*⁻¹ does not exist':
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower == 'dx⁻¹ does not exist':
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower == 'dx^-1 does not exist':
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower == 'dx⁻¹ does not exist':
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower == 'dx⁻¹ does not exist':
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower == 'dx⁻¹ does not exist':
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower == 'dx⁻¹ does not exist':
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower.endswith('does not exist') and 'dx' in lower:
        return rf'\noindent The inverse {inline_math("dX^{-1}")} does not exist in this case.\par'

    if lower.startswith('in si units'):
        return rf'\noindent\emph{{Note.}} {ensure_sentence(escape_tex_text(s))}\par'

    if lower.startswith('put a minus sign in front of each'):
        return r'\noindent Conjugation changes the sign of each basis vector while preserving scalar coefficients.\par'

    if lower.startswith('reverse multiplied order of each'):
        return r'\noindent Reversion reverses the order of each basis-vector product.\par'

    if lower == 'conjugate operation ( )*':
        return r'\noindent\textit{The conjugation signs are summarized below.}\par'

    if lower == 'reverse operation ( )ᴿ':
        return r'\noindent\textit{The reversion signs are summarized below.}\par'

    if lower == 'conjugate reverse operation ( )*ᴿ':
        return r'\noindent\textit{The combined conjugation--reversion signs are summarized below.}\par'

    if lower.startswith('split any vector r into components'):
        return (
            r'\noindent Split the vector '
            + inline_math(r'\mathbf r')
            + r' into parts parallel and perpendicular to the boost direction.\par'
        )

    if lower == 'parallel and perpendicular to u.':
        return (
            r'\noindent The decomposition is taken relative to the unit direction '
            + inline_math(r'\mathbf u')
            + r'.\par'
        )

    if lower.startswith('for n and m in'):
        return (
            r'\noindent For '
            + inline_math(r'n,m \in \{1,2,3\}')
            + r', the following identities hold.\par'
        )

    interval_match = re.match(r'^if (.+), the interval is "(time-like|light-like|space-like)"$', lower)
    if interval_match:
        cond = interval_match.group(1).strip()
        label = interval_match.group(2)
        return rf'\noindent If {inline_math(cond)}, then the interval is \emph{{{label}}}.\par'

    if lower.startswith('if ') and ', the inverse does not exist' in lower:
        cond = s[3:].split(',', 1)[0].strip()
        return rf'\noindent If {inline_math(cond)}, then the inverse does not exist.\par'

    if lower.startswith('if ') and inline_rhs_looks_math(s[3:]):
        cond = s[3:]
        if ' and ' in cond:
            parts = [part.strip() for part in cond.split(' and ')]
            rendered = ' and '.join(inline_math(part) for part in parts)
            return rf'\noindent If {rendered}, then the following relation holds.\par'
        return rf'\noindent If {inline_math(cond)}, then the following relation holds.\par'

    if lower.startswith('for ') and inline_rhs_looks_math(s[4:]):
        return rf'\noindent For {inline_math(s[4:])}, the following formulas apply.\par'

    if lower.startswith('with ') and inline_rhs_looks_math(s[5:]):
        return rf'\noindent With {inline_math(s[5:])}, the formulas reduce to\par'

    if lower.startswith('solve ') and ' in first equation' in lower:
        target = s[6:].split(' in first equation', 1)[0].strip()
        return (
            r'\noindent Solving the first equation for '
            + inline_math(target)
            + r' gives the following relation.\par'
        )

    if lower.startswith('solving ') and ' for ' in lower:
        body = s.rstrip(' ,.')
        parts = body.split(' for ', 1)
        if len(parts) == 2:
            lhs, rhs = parts
            return (
                r'\noindent Solving '
                + inline_math(lhs[len('solving '):] if lhs.lower().startswith('solving ') else lhs)
                + r' for '
                + inline_math(rhs)
                + r' gives the following relation.\par'
            )

    if lower.startswith('components of '):
        expr = s[len('components of '):].strip()
        return rf'\noindent\textit{{Components of {inline_math(expr)}:}}\par'

    label_match = re.match(r'^(real scalar|img scalar|real vector|img vector|scalar|vector|bivector|bi-vector|pseudo-scalar|pseudoscalar)\s*:\s*(.+)$', s, re.IGNORECASE)
    if label_match:
        label = label_match.group(1).strip().capitalize()
        rhs = label_match.group(2).strip()
        if re.fullmatch(r'[A-Za-z \-]+', rhs):
            return rf'\noindent\textit{{{escape_tex_text(label)}:}} {ensure_sentence(escape_tex_text(rhs))}\par'
        return rf'\noindent\textit{{{escape_tex_text(label)}:}} {inline_math(rhs)}.\par'

    decl_match = re.match(r'^(.+?)\s*=\s*([A-Za-z0-9 \-"]+)$', s)
    if decl_match:
        lhs = decl_match.group(1).strip()
        desc = decl_match.group(2).strip()
        if (
            re.fullmatch(r'[A-Za-z0-9 \-"]+', desc)
            and re.search(r'[a-z]', desc)
            and len(desc) > 2
            and ' if ' not in lhs.lower()
            and not inline_rhs_looks_math(desc)
        ):
            return rf'\noindent Let {inline_math(lhs)} denote {escape_tex_text(declaration_phrase(desc))}.\par'

    if re.fullmatch(r'[A-Za-z0-9 ,.\-–—&\'"()]+', s):
        sentence = ensure_sentence(s[0].upper() + s[1:] if s else s)
        return rf'\noindent\textit{{{escape_tex_text(sentence)}}}\par'

    return None


def split_for_alignment(eq: str) -> str:
    s = eq.strip()
    if s.startswith('='):
        return '&' + s
    if s.startswith('+') or s.startswith('-'):
        return r'&\quad ' + s
    for op, aligned in ((':=', '&:='), ('=', '&='), (r'\Rightarrow', r'&\Rightarrow')):
        idx = s.find(op)
        if idx > 0:
            lhs = s[:idx].rstrip()
            rhs = s[idx + len(op):].lstrip()
            return f'{lhs} {aligned} {rhs}'
    return s


def is_continuation_eq(eq: str) -> bool:
    s = eq.strip()
    return s.startswith('=') or s.startswith('+') or s.startswith('-')


def punctuate_equation(eq: str, next_eq: str | None) -> str:
    s = eq.rstrip()
    if re.search(r'[.,;:]$', s):
        return s
    if next_eq is not None and is_continuation_eq(next_eq):
        return s
    if next_eq is None:
        return s + '.'
    return s + ','


def canonical_component_label(label: str) -> str:
    normalized = ' '.join(label.strip().lower().split())
    return normalized.replace('bi-vector', 'bivector').replace('pseudoscalar', 'pseudo-scalar')


def display_component_label(label: str) -> str:
    return COMPONENT_DISPLAY_LABELS.get(label, label.title())


def component_label_tex(label: str) -> str:
    return escape_tex_text(display_component_label(label))


def normalize_component_row_expr(label: str, expr: str) -> str:
    compact = re.sub(r'\s+', '', expr)
    if label.endswith('vector') and compact == '0':
        return r'\mathbf 0'
    return expr


def split_component_annotation(text: str) -> tuple[str, str | None]:
    for marker in ('←', '<-', '->'):
        if marker in text:
            expr, note = text.split(marker, 1)
            return expr.strip(), note.strip()
    return text.strip(), None


def format_component_annotation(note: str) -> str:
    cleaned = ' '.join(note.strip().split())
    if not cleaned:
        return ''

    lower = cleaned.lower()
    if lower.startswith('different than '):
        target = cleaned[len('different than '):].strip()
        return r'\(\leftarrow\)\ \eqtext{different from } \(' + convert_expr(target) + r'\)'
    if lower.startswith('different from '):
        target = cleaned[len('different from '):].strip()
        return r'\(\leftarrow\)\ \eqtext{different from } \(' + convert_expr(target) + r'\)'
    return r'\(\leftarrow\)\ \eqtext{' + escape_tex_text(cleaned) + '}'


def is_zero_component_expr(parts: list[str]) -> bool:
    text = ''.join(parts)
    compact = re.sub(r'\s+', '', text)
    return compact in {'0', '(0)'}


def consume_component_block(items: list[str], start: int) -> tuple[list[str], int] | None:
    stripped = items[start].strip()
    lower = stripped.lower()
    if not lower.startswith('components of '):
        return None

    expr = stripped[len('components of '):].strip()
    rows: list[dict[str, list[str] | str]] = []
    pending_label: str | None = None
    i = start + 1

    def append_to_last(piece: str) -> bool:
        if not rows:
            return False
        expr_piece, note = split_component_annotation(piece)
        if expr_piece:
            rows[-1]['expr_parts'].append(expr_piece)
        if note:
            rows[-1]['notes'].append(note)
        return True

    while i < len(items):
        raw = items[i]
        candidate = raw.strip()
        if not candidate:
            break
        if candidate.lower().startswith('components of '):
            break

        label_match = COMPONENT_ROW_RE.match(candidate)
        if label_match:
            pending_label = canonical_component_label(label_match.group(1))
            rhs = label_match.group(2).strip()
            if rhs:
                expr_piece, note = split_component_annotation(rhs)
                rows.append({
                    'label': pending_label,
                    'expr_parts': [expr_piece] if expr_piece else [],
                    'notes': [note] if note else [],
                })
                pending_label = None
            i += 1
            continue

        if pending_label is not None:
            expr_piece, note = split_component_annotation(candidate)
            rows.append({
                'label': pending_label,
                'expr_parts': [expr_piece] if expr_piece else [],
                'notes': [note] if note else [],
            })
            pending_label = None
            i += 1
            continue

        if candidate.startswith(('←', '<-', '->')):
            if append_to_last(candidate):
                i += 1
                continue
            break

        if raw[: len(raw) - len(raw.lstrip())]:
            if append_to_last(candidate):
                i += 1
                continue

        break

    if not rows:
        return None

    global COMPONENT_BLOCK_COUNTER
    COMPONENT_BLOCK_COUNTER += 1
    block_label = f'tab:components-{COMPONENT_BLOCK_COUNTER:04d}'
    rendered: list[str] = [rf'\begin{{componenttable}}[{block_label}]{{{convert_expr(expr)}}}']
    row_map = {str(row['label']): row for row in rows}
    split_mode = any(label in {'real scalar', 'img scalar', 'real vector', 'img vector'} for label in row_map)

    ordered_labels: list[str]
    if split_mode:
        ordered_labels = ['real scalar', 'img scalar', 'real vector', 'img vector']
    else:
        ordered_labels = [label for label in ('scalar', 'vector') if label in row_map]
        if not ordered_labels:
            ordered_labels = list(row_map.keys())

    generic_scalar_zero = 'scalar' in row_map and is_zero_component_expr(row_map['scalar']['expr_parts'])
    generic_vector_zero = 'vector' in row_map and is_zero_component_expr(row_map['vector']['expr_parts'])

    for label in ordered_labels:
        row = row_map.get(label)
        if row is None and split_mode:
            if label.endswith('scalar') and generic_scalar_zero:
                row = {'label': label, 'expr_parts': ['0'], 'notes': []}
            elif label.endswith('vector') and generic_vector_zero:
                row = {'label': label, 'expr_parts': ['0'], 'notes': []}
            else:
                row = {'label': label, 'expr_parts': ['0'], 'notes': []}

        if row is None:
            continue

        expr_parts = [convert_expr(part) for part in row['expr_parts'] if part]
        row_expr = ' '.join(expr_parts) if expr_parts else '0'
        row_expr = re.sub(r'[.,;:]\s*$', '', row_expr.strip())
        row_expr = normalize_component_row_expr(label, row_expr)
        row_notes = ' '.join(
            formatted for formatted in (format_component_annotation(note) for note in row['notes']) if formatted
        )
        rendered.append(
            rf'\componentrow{{{component_label_tex(label)}}}{{{row_expr}}}{{{row_notes}}}'
        )
    rendered.append(r'\end{componenttable}')
    return rendered, i


def consume_recovery_table(items: list[str], start: int) -> tuple[list[str], int] | None:
    stripped = items[start].strip()
    if stripped.lower() != 'recovering elements':
        return None

    rows: list[dict[str, list[str] | str]] = []
    pending_label: str | None = None
    i = start + 1

    def append_to_last(piece: str) -> bool:
        if not rows:
            return False
        expr_piece, note = split_component_annotation(piece)
        if expr_piece:
            rows[-1]['expr_parts'].append(expr_piece)
        if note:
            rows[-1]['notes'].append(note)
        return True

    while i < len(items):
        raw = items[i]
        candidate = raw.strip()
        if not candidate:
            break
        if format_prose_line(candidate) is not None and not COMPONENT_ROW_RE.match(candidate):
            break

        label_match = COMPONENT_ROW_RE.match(candidate)
        if label_match:
            pending_label = canonical_component_label(label_match.group(1))
            rhs = label_match.group(2).strip()
            if rhs:
                expr_piece, note = split_component_annotation(rhs)
                rows.append({
                    'label': pending_label,
                    'expr_parts': [expr_piece] if expr_piece else [],
                    'notes': [note] if note else [],
                })
                pending_label = None
            i += 1
            continue

        if pending_label is not None:
            expr_piece, note = split_component_annotation(candidate)
            rows.append({
                'label': pending_label,
                'expr_parts': [expr_piece] if expr_piece else [],
                'notes': [note] if note else [],
            })
            pending_label = None
            i += 1
            continue

        if candidate.startswith(('←', '<-', '->')):
            if append_to_last(candidate):
                i += 1
                continue
            break

        if raw[: len(raw) - len(raw.lstrip())]:
            if append_to_last(candidate):
                i += 1
                continue

        break

    if not rows:
        return None

    global COMPONENT_BLOCK_COUNTER
    COMPONENT_BLOCK_COUNTER += 1
    block_label = f'tab:components-{COMPONENT_BLOCK_COUNTER:04d}'
    rendered = [
        r'\noindent The recovery formulas for the scalar, vector, real, and imaginary parts of the complex paravector are collected in Table~\ref{'
        + block_label
        + r'}.\par',
        rf'\begin{{componenttable}}[{block_label}]{{Z}}',
    ]
    ordered_labels = [
        'scalar',
        'vector',
        'real',
        'img',
        'real scalar',
        'img scalar',
        'real vector',
        'img vector',
    ]
    row_map = {str(row['label']): row for row in rows}
    for label in ordered_labels:
        row = row_map.get(label)
        if row is None:
            continue
        expr_parts = [convert_expr(part) for part in row['expr_parts'] if part]
        row_expr = ' '.join(expr_parts) if expr_parts else '0'
        row_expr = re.sub(r'[.,;:]\s*$', '', row_expr.strip())
        row_expr = normalize_component_row_expr(label, row_expr)
        row_notes = ' '.join(
            formatted for formatted in (format_component_annotation(note) for note in row['notes']) if formatted
        )
        rendered.append(
            rf'\componentrow{{{component_label_tex(label)}}}{{{row_expr}}}{{{row_notes}}}'
        )
    rendered.append(r'\end{componenttable}')
    return rendered, i


def consume_grade_component_table(items: list[str], start: int) -> tuple[list[str], int] | None:
    stripped = items[start].strip()
    lower = stripped.lower().rstrip('.')
    if lower not in {
        'every product can be reduced to the following grade components',
        'all product components can be reduced to',
    }:
        return None

    rows: list[dict[str, str]] = []
    i = start + 1
    while i < len(items):
        candidate = items[i].strip()
        if not candidate:
            break
        label_match = COMPONENT_ROW_RE.match(candidate)
        if not label_match:
            break
        label = canonical_component_label(label_match.group(1))
        if label not in GRADE_LABELS:
            break
        rhs = label_match.group(2).strip()
        rows.append({'label': label, 'expr': rhs})
        i += 1

    if not rows:
        return None

    global GRADE_TABLE_COUNTER
    GRADE_TABLE_COUNTER += 1
    block_label = f'tab:sigma-grade-components-{GRADE_TABLE_COUNTER:02d}'
    rendered = [
        r'\noindent The basic grade decomposition of the sigma algebra is summarized in Table~\ref{'
        + block_label
        + r'}, which is the classification used throughout the later product expansions.\par',
        rf'\begin{{gradetable}}[{block_label}]{{Grade decomposition of the sigma algebra.}}',
    ]
    ordered_labels = ['scalar', 'vector', 'bivector', 'pseudo-scalar']
    row_map = {row['label']: row for row in rows}
    for label in ordered_labels:
        row = row_map.get(label)
        if row is None:
            continue
        rendered.append(
            rf'\graderow{{{escape_tex_text(display_component_label(label))}}}{{{convert_expr(row["expr"])}}}'
        )
    rendered.append(r'\end{gradetable}')
    return rendered, i


def is_auxiliary_chain_line(eq: str) -> bool:
    s = eq.strip()
    return s.startswith(r'\qquad') or s.startswith(r'\quad')


def partition_math_lines(lines: list[str]) -> list[list[str]]:
    if not lines:
        return []

    groups: list[list[str]] = [[lines[0]]]
    saw_continuation = False

    for eq in lines[1:]:
        if is_continuation_eq(eq):
            groups[-1].append(eq)
            saw_continuation = True
            continue

        if is_auxiliary_chain_line(eq):
            groups[-1].append(eq)
            continue

        if saw_continuation:
            groups.append([eq])
            saw_continuation = False
            continue

        groups[-1].append(eq)

    return groups


def math_group_kind(lines: list[str]) -> str:
    if len(lines) == 1:
        return 'single'
    if any(is_continuation_eq(eq) for eq in lines[1:]):
        return 'derivation'
    return 'identity_list'


def render_math_group(lines: list[str], has_recent_prose: bool) -> list[str]:
    if not lines:
        return []

    kind = math_group_kind(lines)
    out: list[str] = []

    if kind == 'single':
        if not has_recent_prose:
            out.append(r'\noindent The following numbered relation is recorded for later use.\par')
        out.append(r'\begin{equation}')
        out.append(punctuate_equation(lines[0], None))
        out.append(r'\end{equation}')
        return out

    if kind == 'derivation':
        out.append(r'\noindent The intermediate steps are written out explicitly as follows.\par')
        out.append(r'\begin{align}')
        for idx, eq in enumerate(lines):
            next_eq = lines[idx + 1] if idx < len(lines) - 1 else None
            rendered = punctuate_equation(split_for_alignment(eq), split_for_alignment(next_eq) if next_eq else None)
            if idx < len(lines) - 1:
                rendered += r'\notag \\'
            out.append(rendered)
        out.append(r'\end{align}')
        out.append(r'\noindent The final numbered line is the derived relation used in the subsequent discussion.\par')
        return out

    if not has_recent_prose:
        out.append(r'\noindent The following numbered identities fix the notation and formulas used in this block.\par')
    out.append(r'\begin{align}')
    for idx, eq in enumerate(lines):
        next_eq = lines[idx + 1] if idx < len(lines) - 1 else None
        rendered = punctuate_equation(split_for_alignment(eq), split_for_alignment(next_eq) if next_eq else None)
        suffix = r' \\' if idx < len(lines) - 1 else ''
        out.append(rendered + suffix)
    out.append(r'\end{align}')
    return out


def render_item_lines(items: list[str]) -> list[str]:
    out: list[str] = []
    i = 0
    has_recent_prose = False
    while i < len(items):
        raw = items[i]
        stripped = raw.strip()
        if not stripped:
            out.append('\\medskip')
            has_recent_prose = False
            i += 1
            continue

        recovery_table = consume_recovery_table(items, i)
        if recovery_table is not None:
            rendered, next_index = recovery_table
            out.extend(rendered)
            has_recent_prose = False
            i = next_index
            continue

        grade_table = consume_grade_component_table(items, i)
        if grade_table is not None:
            rendered, next_index = grade_table
            out.extend(rendered)
            has_recent_prose = False
            i = next_index
            continue

        component_block = consume_component_block(items, i)
        if component_block is not None:
            rendered, next_index = component_block
            out.extend(rendered)
            has_recent_prose = False
            i = next_index
            continue

        prose = format_prose_line(stripped)
        if prose is not None:
            out.append(prose)
            note = qualification_note(stripped)
            if note:
                out.append(rf'\noindent\flag{{{note}}}\par')
            has_recent_prose = True
            i += 1
            continue

        if is_display_math_line(stripped):
            math_lines: list[str] = []
            notes: list[str] = []
            while i < len(items):
                candidate = items[i].strip()
                if (
                    not candidate
                    or format_prose_line(candidate) is not None
                    or not is_display_math_line(candidate)
                ):
                    break
                math_lines.append(convert_expr(candidate))
                note = qualification_note(candidate)
                if note:
                    notes.append(note)
                i += 1
            for group in partition_math_lines(math_lines):
                out.extend(render_math_group(group, has_recent_prose))
                has_recent_prose = False
            for note in dict.fromkeys(notes):
                out.append(rf'\noindent\flag{{{note}}}\par')
            continue

        text = escape_tex_text(stripped)
        indent = indent_tex(raw)
        out.append(f'\\noindent {indent}{text}\\par')
        note = qualification_note(stripped)
        if note:
            out.append(rf'\noindent\flag{{{note}}}\par')
        has_recent_prose = True
        i += 1

    return out


def emit_equation_catalog(blocks: list[tuple[str, list[tuple[int, str]]]]) -> tuple[list[str], int, int]:
    total = sum(len(items) for _, items in blocks)

    out = []
    out.append('% Auto-generated equation catalog. Do not hand-edit.')
    out.append('\\section{Equation Line Catalog}')
    out.append('This generated catalog converts every equation-bearing line from \\texttt{README.txt} into LaTeX math, preserving source order and line numbers.')
    out.append(f'\\textbf{{Total converted equation lines: {total}.}}')

    block_no = 0
    ref_no = 0
    for title, items in blocks:
        if not items:
            continue
        block_no += 1
        title_tex = escape_tex_text(title)
        out.append('')
        out.append(f'\\subsection*{{Block {block_no:03d}: {title_tex}}}')
        for line_no, src in items:
            ref_no += 1
            eq = convert_expr(src)
            # Keep original source line visible for auditability.
            src_tex = escape_tex_text(src)
            out.append(f'\\paragraph{{R-{ref_no:04d} (README line {line_no})}}')
            out.append(f'\\phantomsection\\label{{eq:readme-r{ref_no:04d}}}')
            out.append(f'\\textit{{Source:}} \\texttt{{{src_tex}}}')
            out.append('\\[')
            out.append(eq)
            out.append('\\]')

    return out, total, block_no


def emit_full_conversion(lines: list[str]) -> list[str]:
    out = []
    out.append('% Auto-generated source conversion. Do not hand-edit.')
    out.append('\\section{Integrated Full README Conversion}\\label{sec:full-readme}')
    out.append('This section incorporates every line of \\texttt{README.txt} into the main body. Equation-bearing lines are converted to LaTeX math; prose lines are rendered as structured text; blank and separator lines are explicitly labeled so the full source is accounted for.')
    out.append(f'\\textbf{{Total source lines integrated: {len(lines)}.}}')
    out.append('\\begingroup')
    out.append('\\small')

    expecting_title = False

    for idx, raw in enumerate(lines, start=1):
        line = raw.rstrip('\n')
        stripped = line.strip()

        if re.fullmatch(r'_+', stripped):
            out.append('')
            out.append(rf'\noindent\textbf{{README line {idx}.}} \hrulefill\par')
            expecting_title = True
            continue

        if expecting_title and stripped:
            title_tex = escape_tex_text(stripped)
            out.append(f'\\subsection*{{README line {idx}: {title_tex}}}')
            expecting_title = False
            continue

        if not stripped:
            out.append(rf'\noindent\textbf{{README line {idx}.}} \textit{{[blank line]}}\par')
            continue

        if is_equation_line(stripped):
            eq = convert_expr(stripped)
            out.append(f'\\noindent\\textbf{{README line {idx}.}}')
            out.append('\\[')
            out.append(eq)
            out.append('\\]')
            continue

        text = escape_tex_text(stripped)
        indent = indent_tex(line)
        out.append(f'\\noindent\\textbf{{README line {idx}.}} {indent}{text}\\par')

    out.append('\\endgroup')
    return out


def parse_full_blocks(lines: list[str]) -> list[tuple[str, list[str]]]:
    blocks: list[tuple[str, list[str]]] = []
    current_title = 'Initial Definitions'
    current_items: list[str] = []
    expecting_title = False

    def flush() -> None:
        nonlocal current_items, current_title
        if current_items:
            blocks.append((current_title, current_items))
            current_items = []

    for raw in lines:
        line = raw.rstrip('\n')
        stripped = line.strip()

        if re.fullmatch(r'_+', stripped):
            flush()
            expecting_title = True
            continue

        if expecting_title and stripped:
            current_title = stripped
            expecting_title = False
            continue

        current_items.append(line)

    flush()
    return blocks


def emit_extended_appendix(blocks: list[tuple[str, list[str]]]) -> list[str]:
    out = []
    out.append('% Auto-generated appendix. Do not hand-edit.')
    out.append('\\section{Extended Algebraic and Physical Notes}\\label{app:extended-notes}')
    out.append('This appendix collects the broader catalog of identities, definitions, product expansions, and auxiliary remarks supporting the main text. Material is grouped by topic and rendered in manuscript notation so the full development remains available without interrupting the main exposition. Repeated developments are retained when they support a distinct derivational context, and statements known to require qualification are marked in bold red at the point where they appear.')
    out.append('\\begingroup')
    out.append('\\small')

    for title, items in blocks:
        filtered = [line for line in items if line.strip()]
        if not filtered:
            continue

        title_tex = escape_tex_text(title)
        out.append('')
        out.append(f'\\subsection*{{{title_tex}}}')
        out.extend(render_item_lines(items))

    out.append('\\endgroup')
    return out


def emit_integrated_body(blocks: list[tuple[str, list[str]]]) -> tuple[list[str], int]:
    remaining = list(blocks)
    for skip_title in SKIP_INTEGRATED_TITLES:
        skipped = consume_block(remaining, skip_title)
        while skipped is not None:
            skipped = consume_block(remaining, skip_title)

    out = []
    out.append('% Auto-generated integrated development. Do not hand-edit.')
    out.append('\\section{Integrated Algebraic and Physical Development}\\label{sec:integrated-development}')
    out.append(
        'This section absorbs the broader algebraic and physical catalog into the main body of the manuscript. '
        'The opening summary page and foundational definitions already cover the initial overview material, so the development below begins with the remaining detailed notes and arranges them thematically after the prerequisite concepts have been introduced in the main exposition.'
    )
    out.append(
        'Each subsection therefore serves a double purpose: it preserves the fuller equation inventory while also giving the reader a coherent path from the compact formulas of the earlier sections to the more exhaustive identities, intermediate expansions, and contextual remarks that support them.'
    )

    integrated_count = 0
    for group in INTEGRATED_GROUPS:
        out.append('')
        out.append(f"\\subsection{{{group['title']}}}\\label{{{group['label']}}}")
        out.append(group['intro'])
        for title in group['blocks']:
            block = consume_block(remaining, title)
            if block is None:
                continue
            integrated_count += 1
            raw_title, items = block
            display_title = escape_tex_text(integrated_display_title(raw_title, group['title']))
            intro = integrated_block_intro(raw_title, group['title'])
            out.append('')
            out.append(f'\\subsubsection*{{{display_title}}}')
            out.append(intro)
            out.extend(render_item_lines(items))

    if remaining:
        out.append('')
        out.append('\\subsection{Additional Consolidated Notes}')
        out.append(
            'The following blocks were not assigned to one of the main thematic subsections above, but they are retained here so that the integrated development remains complete.'
        )
        for raw_title, items in remaining:
            integrated_count += 1
            display_title = escape_tex_text(raw_title)
            out.append('')
            out.append(f'\\subsubsection*{{{display_title}}}')
            out.append(
                'This subsection preserves a remaining block whose content is still relevant to the integrated development.'
            )
            out.extend(render_item_lines(items))

    return out, integrated_count


def emit_integrated_group_content(
    group_title: str,
    block_specs: list[str],
    remaining: list[tuple[str, list[str]]],
) -> tuple[list[str], int]:
    out: list[str] = []
    integrated_count = 0

    for title in block_specs:
        block = consume_block(remaining, title)
        if block is None:
            continue
        integrated_count += 1
        raw_title, items = block
        display_title = escape_tex_text(integrated_display_title(raw_title, group_title))
        intro = integrated_block_intro(raw_title, group_title)
        out.append('')
        out.append(f'\\subsubsection*{{{display_title}}}')
        out.append(intro)
        out.extend(render_item_lines(items))

    return out, integrated_count


def emit_integrated_group_files(blocks: list[tuple[str, list[str]]]) -> tuple[dict[str, list[str]], int]:
    remaining = list(blocks)
    for skip_title in SKIP_INTEGRATED_TITLES:
        skipped = consume_block(remaining, skip_title)
        while skipped is not None:
            skipped = consume_block(remaining, skip_title)

    outputs: dict[str, list[str]] = {}
    integrated_count = 0

    for group in INTEGRATED_GROUPS:
        out = ['% Auto-generated integrated section fragment. Do not hand-edit.']
        group_lines, count = emit_integrated_group_content(group['title'], group['blocks'], remaining)
        out.extend(group_lines)
        outputs[group['slug']] = out
        integrated_count += count

    if remaining:
        out = ['% Auto-generated integrated section fragment. Do not hand-edit.']
        out.append('\\subsubsection*{Additional Consolidated Notes}')
        out.append(
            'The following blocks were not assigned to one of the primary section-local groupings above, but they are retained here so the integrated manuscript remains complete.'
        )
        for raw_title, items in remaining:
            integrated_count += 1
            display_title = escape_tex_text(raw_title)
            out.append('')
            out.append(f'\\paragraph{{{display_title}}}')
            out.append(
                'This remaining block is preserved here because it still contributes to the broader representation and extension material.'
            )
            out.extend(render_item_lines(items))
        outputs['representations'].extend(out[1:])

    return outputs, integrated_count


def finalize_generated_tex(lines: list[str]) -> list[str]:
    text = '\n'.join(lines)
    text = text.replace('\u200b', '')
    text = text.replace('₊', '_{+}')
    text = text.replace('₋', '_{-}')
    text = text.replace('2×2', r'$2\times2$')
    text = text.replace('electro-magnetic', 'electromagnetic')
    text = text.replace(r'\\in ', r'\in ')
    text = text.replace(r'\mathbb{\mathbf C}^{3}', r'\mathbb{C}^{3}')
    text = text.replace(r'\mathbb{\mathbf C}', r'\mathbb{C}')
    text = text.replace(r'\partial _{\mathbf r}', r'\delr{}')
    text = text.replace(r'\partial_{\mathbf r}', r'\delr{}')
    text = text.replace(r'\partial _{r}', r'\delr{}')
    text = text.replace(r'\partial_{r}', r'\delr{}')
    text = text.replace(r'\partial _{\mathbf v}', r'\delv{}')
    text = text.replace(r'\partial_{\mathbf v}', r'\delv{}')
    text = text.replace(r'\partial _{v}', r'\delv{}')
    text = text.replace(r'\partial_{v}', r'\delv{}')
    text = text.replace(r'\nabla', r'\gradv')
    text = text.replace(r'\mathbb{R}+\mathbb{R}^{3}', r'\mathbb{R}\oplus\mathbb{R}^{3}')
    text = text.replace(r'\mathbb{C}+\mathbb{C}^{3}', r'\mathbb{C}\oplus\mathbb{C}^{3}')
    text = text.replace(r'f(+1)', r'f(+)')
    text = text.replace(r'f(-1)', r'f(-)')
    text = text.replace(r'{1,2,3}, {2,3,1}, {3,1,2}', r'\{1,2,3\}, \{2,3,1\}, \{3,1,2\}')
    text = text.replace(r'{3,2,1}, {1,3,2}, {2,1,3}', r'\{3,2,1\}, \{1,3,2\}, \{2,1,3\}')
    text = text.replace(
        r'\sigv = \left\{\sigma _{1}, \sigma _{2}, \sigma _{3}\right\}\qquad \text{symbolic vector}.',
        r'\sigv = \{ \sigma _{1}, \sigma _{2}, \sigma _{3} \}\qquad \eqtext{vector of algebraic objects}.',
    )
    text = text.replace(
        r'\sigv = (\sigma _{1}, \sigma _{2}, \sigma _{3})\qquad \text{vector of algebraic objects}.',
        r'\sigv = \{ \sigma _{1}, \sigma _{2}, \sigma _{3} \}\qquad \eqtext{vector of algebraic objects}.',
    )
    text = text.replace(r'\text{complex scalar}', r'\eqtext{complex scalar}')
    text = text.replace(r'\text{real scalar}', r'\eqtext{real scalar}')
    text = text.replace(r'\text{complex vector}', r'\eqtext{complex vector}')
    text = text.replace(r'\text{real vector}', r'\eqtext{real vector}')
    text = text.replace(r'\text{vector of algebraic objects}', r'\eqtext{vector of algebraic objects}')
    text = text.replace(r'\text{real scalar + complex vector}', r'\eqtext{real scalar + complex vector}')
    text = text.replace(r'\text{real scalar + real vector}', r'\eqtext{real scalar + real vector}')
    text = text.replace(r'\qquad different than ', r'\qquad \text{different from }')
    text = text.replace(r'\qquad 4 force.', r'\qquad \text{4-force}.')
    text = text.replace(r'\qquad 4 force,', r'\qquad \text{4-force},')
    text = text.replace(r'&= 0 else.', r'&= 0 \text{ else}.')
    text = text.replace(r'1 &= vector,', r'1 &= \text{vector},')
    text = text.replace(r'2 &= bi-vector,', r'2 &= \text{bi-vector},')
    text = text.replace(r'3 &= pseudo-scalar,', r'3 &= \text{pseudo-scalar},')
    text = text.replace(r'action &=', r'\text{action} &=')
    text = text.replace(r'\left\{', r'\{')
    text = text.replace(r'\right\}', r'\}')
    text = re.sub(r'\\Sigma\s*_\{([A-Za-z0-9]+)\}', r'\\sum_{\1}', text)
    text = re.sub(r'\\sqrt\{\}\s*(\\square)', r'\\sqrt{\1}', text)
    text = re.sub(r'\\sqrt\{\}\s*(\\lambda\s*_\{[^}]+\})', r'\\sqrt{\1}', text)
    text = re.sub(r'\\sqrt\{\}\s*(\\mathbf\s+[A-Za-z](?:\^\{[^}]+\})?)', r'\\sqrt{\1}', text)
    text = re.sub(r'\\sqrt\{\}\s*([A-Za-z](?:_\{[^}]+\})?(?:\^\{[^}]+\})?)', r'\\sqrt{\1}', text)
    text = text.replace(
        r'\langle \psi | \psi \rangle = expectation of \psi ^{\mathrm{R}} \psi = 1.',
        r'\langle \psi | \psi \rangle = \operatorname{exp}(\psi^{\mathrm{R}}\psi) = 1.',
    )
    text = text.replace(
        r'\langle \psi | M | \psi \rangle = expectation of \psi ^{\mathrm{R}} M \psi.',
        r'\langle \psi | M | \psi \rangle = \operatorname{exp}(\psi^{\mathrm{R}} M \psi).',
    )
    text = text.replace(r'\operatorname{Exp}', r'\operatorname{exp}')
    text = text.replace(
        r'\noindent\emph{Note.} ( )ᴿ is equivalent to the Hermetian operator.\par',
        r'\noindent\emph{Note.} In this block representation, \(^{\mathrm R}\) plays the role of Hermitian conjugation.\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} ( )ᴿ is the Hermitian operator.\par',
        r'\noindent\emph{Note.} In this block representation, \(^{\mathrm R}\) plays the role of Hermitian conjugation.\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} E + P • σ = i ∂*.\par',
        r'\noindent\emph{Note.} Equivalently, \(E+\mathbf P\cdot\sigv=i\partial^{*}\).\par',
    )
    # The plain-text README format has no bold Greek support, so collapse the
    # block-field symbol back to ordinary capital Psi during conversion.
    text = text.replace(r'\bPsi', r'\Psi')
    text = text.replace(
        'The formulas below make products of real scalars, vectors, and paravectors explicit and collect identities used later in the manuscript.',
        'This subsection expands the basic scalar, vector, and paravector product rules used throughout the manuscript.',
    )
    text = text.replace(
        r'\noindent\textit{Consider the 2-element vector.}\par',
        r'\noindent\textit{Introduce the two-component block field.}\par',
    )
    text = text.replace(
        r'\Psi = \left\{ \Psi _{1} ; \Psi _{2} \right\}.',
        r'\Psi = \{ \Psi _{1}; \Psi _{2} \}.',
    )
    text = text.replace(
        r'\noindent define the $2\times2$ operator matrix\par',
        r'\noindent\textit{Define the corresponding $2\times2$ block operator.}\par',
    )
    text = text.replace(
        r'W(Z) &:= \left\{ 0 , Z ; Z^{*\mathrm{R}} , 0 \right\}\qquad 2\times 2 matrix, \\',
        r'\Wmat(Z) &:= \{ 0, Z; Z^{*\mathrm{R}}, 0 \}\qquad \text{$2\times2$ block matrix}, \\',
    )
    text = text.replace(
        r'W(Z) &:= \{ 0 , Z ; Z^{*\mathrm{R}} , 0 \}\qquad 2\times 2 matrix, \\',
        r'\Wmat(Z) &:= \{ 0, Z; Z^{*\mathrm{R}}, 0 \}\qquad \text{$2\times2$ block matrix}, \\',
    )
    text = text.replace(
        r'W( i \partial ) &= i { 0 , \partial ; \partial ^{*} , 0 }, \\',
        r'\Wmat( i \partial ) &= i \{ 0, \partial; \partial ^{*}, 0 \}, \\',
    )
    text = text.replace(
        r'\noindent\textit{If a and b are real scalars,.}\par',
        r'\noindent\textit{If \(a,b\in\R\), then}\par',
    )
    text = text.replace(
        r'\noindent\textit{And its mass conjugate.}\par',
        r'\noindent\textit{Introduce also the mass-conjugate operator.}\par',
    )
    text = text.replace(
        "\\noindent using the scalar solution Ψ' from KG\\par",
        r'\noindent\textit{Using a scalar solution \(\psi^\prime\) of the Klein--Gordon equation, one may construct a Dirac solution by}\par',
    )
    text = text.replace(
        r'\noindent With \(a = 1 and b = 0\), the formulas reduce to\par',
        r'\noindent\textit{With \(a=1\) and \(b=0\), the formulas reduce to}\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} definitions differ across sources.\par',
        r'\noindent\emph{Note.} Sign conventions for \(F\) vary across the literature; the present manuscript keeps the convention fixed by Eqs.\ \eqref{eq:F-Phi} and \eqref{eq:maxwell}.\par',
    )
    text = text.replace(
        r'\noindent using paravector metric <<Z>>² = Z Z*ᴿ\par',
        r'\noindent\textit{Using the paravector metric } \(\langle\!\langle Z\rangle\!\rangle^{2}=Z Z^{*\mathrm{R}}\)\textit{, one finds}\par',
    )
    text = text.replace(
        r"\noindent ≥ ( u × u' )²\par",
        r'\noindent\textit{Hence the lower bound is governed by } \((\mathbf u\times \mathbf u^{\prime})^{2}\)\textit{.}\par',
    )
    text = text.replace(r'\mathbf q &= (q_{1}, q_{2}, q_{3}), \\', r'\mathbf q &= \{ q_{1}, q_{2}, q_{3} \}, \\')
    text = text.replace(r'\mathbf b &= (b_{1}, b_{2}, b_{3}), \\', r'\mathbf b &= \{ b_{1}, b_{2}, b_{3} \}, \\')
    text = text.replace(r'\mathbf u &= (u_{1}, u_{2}, u_{3}), \mathbf u^{2} = 1, \\', r'\mathbf u &= \{ u_{1}, u_{2}, u_{3} \}, \mathbf u^{2} = 1, \\')
    text = text.replace(r'\mathbf n &= (n_{1}, n_{2}, n_{3}), \mathbf n^{2} = 1, \\', r'\mathbf n &= \{ n_{1}, n_{2}, n_{3} \}, \mathbf n^{2} = 1, \\')
    text = text.replace(r'\mathbf v &= d\mathbf r/dt = (v_{1}, v_{2}, v_{3}) \in \mathbb{R}^{3}, \\', r'\mathbf v &= d\mathbf r/dt = \{ v_{1}, v_{2}, v_{3} \} \in \mathbb{R}^{3}, \\')
    text = text.replace(r'\mathbf u &= e_{3} = (0, 0, 1), \\', r'\mathbf u &= e_{3} = \{ 0, 0, 1 \}, \\')
    text = text.replace(r'\mathbf P \times \mathbf P &= ( [P_{2}, P_{3}], -[P_{1}, P_{3}], [P_{1}, P_{2}] ), \\', r'\mathbf P \times \mathbf P &= \{ [P_{2}, P_{3}], -[P_{1}, P_{3}], [P_{1}, P_{2}] \}, \\')
    text = text.replace(r'\mathbf A &= \text{3-element vector}, \\', r'\mathbf A &\in \mathbb{R}^{3}, \\')
    text = text.replace(r'\mathbf A &= \text{3-element vector}.', r'\mathbf A \in \mathbb{R}^{3}.')
    text = text.replace(r'\mathbf B &= \text{3-element vector}.', r'\mathbf B \in \mathbb{R}^{3}.')
    text = text.replace(r'\mathbf B \in \mathbb{R}^{3}.', r'\mathbf B &\in \mathbb{R}^{3}.')
    text = text.replace(
        r'\noindent\emph{Note.} C² ≠ ‖C‖².\par',
        r'\noindent\emph{Note.} The quadratic product \(C^{2}\) is not the same object as the norm square \(\lVert \mathbf C\rVert^{2}\).\par',
    )
    text = text.replace(
        """\\[
\\mathbf A \\times \\mathbf B = (A_{2}B_{3} - A_{3}B_{2},
\\]
\\noindent \\hspace*{6.0em}–A₁B₃ + A₃B₁,\\par
\\noindent \\hspace*{7.0em}A₁B₂ – A₂B₁)\\par""",
        r"""\[
\mathbf A \times \mathbf B
=
\{ A_{2}B_{3} - A_{3}B_{2}, -A_{1}B_{3} + A_{3}B_{1}, A_{1}B_{2} - A_{2}B_{1} \}.
\]""",
    )
    text = text.replace(
        """\\begin{equation}
\\mathbf A \\times \\mathbf B = (A_{2}B_{3} - A_{3}B_{2},
\\end{equation}
\\noindent \\hspace*{6.0em}–A₁B₃ + A₃B₁,\\par
\\noindent \\hspace*{7.0em}A₁B₂ – A₂B₁)\\par""",
        r"""\begin{equation}
\mathbf A \times \mathbf B
=
\{ A_{2}B_{3} - A_{3}B_{2}, -A_{1}B_{3} + A_{3}B_{1}, A_{1}B_{2} - A_{2}B_{1} \}.
\end{equation}""",
    )
    text = text.replace(
        """\\noindent\\textit{0 1 2 3.}\\par
\\noindent ( )*      +  –  +  –\\par
\\noindent ( )ᴿ      +  +  –  –\\par
\\noindent ( )*ᴿ     +  –  –  +\\par""",
        r"""\[
\begin{array}{c|cccc}
& 0 & 1 & 2 & 3 \\
( )^{*} & + & - & + & - \\
( )^{\mathrm{R}} & + & + & - & - \\
( )^{*\mathrm{R}} & + & - & - & +
\end{array}
\]""",
    )
    text = text.replace(
        r'\noindent XY ∈ ℝ+ℂ³\par',
        r'\noindent The product \(XY\) generally lies in \(\mathbb{R}\oplus\mathbb{C}^{3}\).\par',
    )
    text = text.replace(
        r'\[' '\n' r'\qquad \text{different from }Z^{*}Z.' '\n' r'\]',
        r'\noindent\emph{Note.} This differs from \(Z^{*}Z\).\par',
    )
    text = text.replace(
        r'\[' '\n' r'\qquad \text{different from }ZZ^{*}.' '\n' r'\]',
        r'\noindent\emph{Note.} This differs from \(ZZ^{*}\).\par',
    )
    text = text.replace(
        r"\noindent \hspace*{3.5em}– A × B' – B × A'\par",
        r'\noindent\textit{with the additional cross-product contribution } \(-\mathbf A\times\mathbf B^{\prime}-\mathbf B\times\mathbf A^{\prime}\)\textit{.}\par',
    )
    text = text.replace(
        r"\noindent \hspace*{4.0em}+ A × A' – B × B'\par",
        r'\noindent\textit{with the additional cross-product contribution } \(\mathbf A\times\mathbf A^{\prime}-\mathbf B\times\mathbf B^{\prime}\)\textit{.}\par',
    )
    text = text.replace(
        r'\noindent Z ∈ ℂ+ℂ³\par',
        r'\noindent Here \(Z\in\mathbb{C}\oplus\mathbb{C}^{3}\).\par',
    )
    text = text.replace(
        r"\noindent \hspace*{1.0em}+ f(b) g(a') Π(–)Π(+)' +f(b) g(b') Π(–)Π(–)'\par",
        r'\noindent\textit{followed by the terms } \(f(b)g(a^{\prime})\Pi(-)\Pi(+)^{\prime}+f(b)g(b^{\prime})\Pi(-)\Pi(-)^{\prime}\)\textit{.}\par',
    )
    text = text.replace(
        r"\noindent \hspace*{1.0em}+ f(b) g(a') Π(–)Π(+)' +f(b) g(b') Π(–)Π(–)'\par",
        r'\noindent\textit{followed by the terms } \(f(b)g(a^{\prime})\Pi(-)\Pi(+)^{\prime}+f(b)g(b^{\prime})\Pi(-)\Pi(-)^{\prime}\)\textit{.}\par',
    )
    text = text.replace(
        r"\noindent \hspace*{5.5em}– (dt² – dr²) – (dt'² – dr'²)\par",
        r'\noindent\textit{That is, one subtracts } \((dt^{2}-d\mathbf r^{2})+(dt^{\prime 2}-d\mathbf r^{\prime 2})\)\textit{ from the combined quadratic form.}\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} in SI units, ds = c dτ.\par',
        r'\noindent\emph{Note.} In SI units one has \(ds=c\,d\tau\).\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} if m = 0, E² = P² ← light.\par',
        r'\noindent\emph{Note.} For \(m=0\), the mass-shell relation reduces to \(E^{2}=\mathbf P^{2}\), corresponding to lightlike motion.\par',
    )
    text = text.replace(
        r'\noindent ≈ m v\par',
        r'\noindent\emph{Thus the leading momentum approximation is } \(m\mathbf v\)\emph{.}\par',
    )
    text = text.replace(
        r'\noindent ≈ ½ m v² + m\par',
        r'\noindent\emph{Thus the leading energy approximation is } \(\tfrac{1}{2}m\mathbf v^{2}+m\)\emph{.}\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} with σ-algebra.\par',
        r'\noindent\emph{Note.} This expansion uses the sigma-algebra conventions fixed earlier in the manuscript.\par',
    )
    text = text.replace(
        r'\noindent some use the  (–,+,+,+) signature\par',
        r'\noindent\emph{Note.} Some authors instead use the metric signature \((- ,+,+,+)\).\par',
    )
    text = text.replace(
        r'\noindent some use □²\par',
        r'\noindent\emph{Note.} Some authors denote the wave operator with a different box-symbol convention.\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} □ λ is a scalar field.\par',
        r'\noindent\emph{Note.} The quantity \(\square\lambda\) is again a scalar field.\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} 𝓔 = electric field (not energy).\par',
        r'\noindent\emph{Note.} Here \(\mathbf E\) denotes the electric field, not the energy operator \(E\).\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} 𝓑 does no work.\par',
        r'\noindent\emph{Note.} The magnetic field contributes no work term in the scalar component.\par',
    )
    text = text.replace(
        r'\noindent\textit{If boost or rotor,.}\par',
        r'\noindent\textit{For boosts and rotors, the inverse simplifies to}\par',
    )
    text = text.replace(
        r'\noindent\textit{Conjugate and reverse.}\par',
        r'\noindent\textit{Applying conjugation and reversion gives}\par',
    )
    text = text.replace(
        r'\noindent\textit{E.g., non-inertial frame.}\par',
        r'\noindent\emph{This interpretation also applies in non-inertial settings.}\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} img( F U ) swaps 𝓔 and 𝓑.\par',
        r'\noindent\emph{Note.} The imaginary part of \(FU\) exchanges the electric and magnetic roles and is therefore not the sector used in the physical force law.\par',
    )
    text = text.replace(
        r'\noindent \hspace*{2.0em}+ q V – q v • A\par',
        r'\noindent\textit{Adding the potential term contributes } \(qV-q\,\mathbf v\cdot\mathbf A\)\textit{.}\par',
    )
    text = text.replace(
        """\\noindent \\hspace*{2.0em}+ (–T⁻¹ dT) (T T⁻¹) X T\\par
\\noindent \\hspace*{2.0em}+ T⁻¹ X (T T⁻¹) dT\\par""",
        r'\noindent\textit{The two intermediate terms are } \(( -T^{-1}dT)(TT^{-1})XT\)\textit{ and } \(T^{-1}X(TT^{-1})dT\)\textit{.}\par',
    )
    text = text.replace(
        r'\noindent \hspace*{2.0em}if  n ≠ m\par',
        r'\noindent For \(n\neq m\), the mixed products transform as follows.\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} L⁻¹ = L* = L*ᴿ.\par',
        r'\noindent\emph{Note.} The boost inverse satisfies \(L^{-1}=L^{*}=L^{*\mathrm{R}}\).\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} R⁻¹ = Rᴿ = R*ᴿ.\par',
        r'\noindent\emph{Note.} The rotor inverse satisfies \(R^{-1}=R^{\mathrm{R}}=R^{*\mathrm{R}}\).\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} Ψ is a complex scalar.\par',
        r'\noindent\emph{Note.} Here \(\psi\) is a complex scalar field.\par',
    )
    text = text.replace(r'+ ···.', r'+ \cdots.')
    text = text.replace(
        r'\noindent\emph{Note.} R(φ + 2π) = -R(φ).\par',
        r'\noindent\emph{Note.} The rotor changes sign under a full \(2\pi\) rotation: \(R(\varphi+2\pi)=-R(\varphi)\).\par',
    )
    text = text.replace(
        r"\mathbf r' \cdot \mathbf q = Q (\mathbf r \cdot \mathbf q) Q^{\mathrm{R}}\qquad \text{see R(θ,u)}.",
        r'\mathbf r^{\prime}\cdot\mathbf q = Q(\mathbf r\cdot\mathbf q)Q^{\mathrm{R}}, \qquad \text{cf. Eq.\ \eqref{eq:rotor}}.',
    )
    text = text.replace(r'\text{spin up + spin down}', r'\text{spin decomposition}')
    text = text.replace(
        r'\noindent \hspace*{2.0em}+ (B • dr)(B* • dr)\par',
        r'\noindent\textit{and the remaining quadratic contribution is } \((\mathbf B\cdot d\mathbf r)(\mathbf B^{*}\cdot d\mathbf r)\)\textit{.}\par',
    )
    text = text.replace(
        r'\noindent\emph{Note.} to remain time-like, v² < 1.\par',
        r'\noindent\emph{Note.} The timelike condition requires \(\mathbf v^{2}<1\).\par',
    )
    text = text.replace(
        r'\noindent \hspace*{1.0em}for n ≠ m\par',
        r'\noindent For \(n\neq m\), the mixed matrix products satisfy\par',
    )
    return text.split('\n')


def main() -> None:
    lines = README.read_text(encoding='utf-8').splitlines()
    blocks = parse_blocks(lines)
    full_blocks = parse_full_blocks(lines)
    equation_out, total, block_no = emit_equation_catalog(blocks)
    full_out = emit_full_conversion(lines)
    appendix_out = emit_extended_appendix(full_blocks)
    integrated_out, integrated_count = emit_integrated_body(full_blocks)
    integrated_group_out, integrated_group_count = emit_integrated_group_files(full_blocks)

    equation_out = finalize_generated_tex(equation_out)
    full_out = finalize_generated_tex(full_out)
    appendix_out = finalize_generated_tex(appendix_out)
    integrated_out = finalize_generated_tex(integrated_out)
    integrated_group_out = {
        slug: finalize_generated_tex(content)
        for slug, content in integrated_group_out.items()
    }

    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    GROUP_OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT.write_text('\n'.join(equation_out) + '\n', encoding='utf-8')
    FULL_OUT.write_text('\n'.join(full_out) + '\n', encoding='utf-8')
    APPENDIX_OUT.write_text('\n'.join(appendix_out) + '\n', encoding='utf-8')
    BODY_OUT.write_text('\n'.join(integrated_out) + '\n', encoding='utf-8')
    for slug, content in integrated_group_out.items():
        out_path = GROUP_OUT_DIR / f'{slug}.tex'
        out_path.write_text('\n'.join(content) + '\n', encoding='utf-8')
    print(f'Wrote {OUT} with {total} converted equation lines in {block_no} blocks.')
    print(f'Wrote {FULL_OUT} with {len(lines)} integrated source lines.')
    print(f'Wrote {APPENDIX_OUT} with {len(full_blocks)} topical blocks.')
    print(f'Wrote {BODY_OUT} with {integrated_count} integrated topical blocks.')
    print(f'Wrote {GROUP_OUT_DIR} with {integrated_group_count} integrated topical blocks across section fragments.')


if __name__ == '__main__':
    main()
