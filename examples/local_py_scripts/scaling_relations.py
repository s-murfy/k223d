"""
scaling_relations.py

A pluggable interface to empirical earthquake magnitude <-> rupture-area
(and length) scaling relationships (Wells & Coppersmith 1994, Strasser et
al. 2010, Leonard 2010, Hanks & Bakun 2002, etc).

Design problems this solves
-----------------------------
1. Different papers subdivide earthquakes differently:
     - Wells & Coppersmith (1994): by *faulting style* (strike-slip,
       reverse, normal, or an "all" category).
     - Strasser et al. (2010): by *subduction-zone position* (interface,
       intraslab) -- faulting style isn't relevant to their regressions.
     - Leonard (2010): by *tectonic regime + style*.
     - Hanks & Bakun (2002/2008): no subdivision at all, and the M-logA
       relation is *bilinear* (a break in slope at a hinge magnitude).
   Rather than forcing one shared category system onto every paper, each
   relation declares its OWN valid categories.

2. Papers publish scaling relations in TWO distinct regression directions,
   and these are NOT algebraically interchangeable:
     - "area-on-magnitude":     log10(A) = a + b * M   (e.g. W&C Table 2A)
     - "magnitude-on-area":     M = a + b * log10(A)   (e.g. W&C Table 2B)
   Ordinary least-squares minimizes residuals in the DEPENDENT variable
   only, so a fit of A on M and a fit of M on A are generally different
   lines with different coefficients and different uncertainty (sigma is
   in log10(area) units for the first form, in Mw units for the second).
   This module stores whichever direction(s) a paper actually published,
   and automatically uses the statistically correct one for whichever
   operation you're doing -- falling back to an algebraic inversion of
   the other direction (with a warning) only if the paper didn't publish
   a fit in the direction you need.

3. Some relations regress against Mw directly; others regress against
   log10(seismic moment, M0) instead -- see the moment-conversion section
   below. Either way, this module's user-facing in/out is always Mw.

Core equations
--------------
    Area-on-magnitude form:   log10(A) = a + b * M
    Magnitude-on-area form:   M = a + b * log10(A)

    where
        A = rupture area, in m^2 (this module's native/default unit --
            see the unit-handling section below; can be requested in
            km^2 instead)
        M = moment magnitude, Mw (dimensionless)
        a = regression intercept (specific to the relation, category,
            AND which of the two forms above)
        b = regression slope     (specific to the relation, category,
            AND which of the two forms above)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import math
import warnings


# ---------------------------------------------------------------------------
# Unit handling
# ---------------------------------------------------------------------------
# Different papers fit their regressions in different area units -- most
# large-earthquake area relations (Wells & Coppersmith, Strasser, Leonard)
# use km^2, but some are quoted in m^2. Silently mixing these up doesn't
# error -- it just gives an answer wrong by a factor of 1e6 (1 km^2 =
# 1e6 m^2). To avoid that, every ScalingCoefficients records the unit it
# was FIT in, and conversion happens once, at the boundary, in compute().
# This module's own default in/out unit is m^2, matching k223d's own
# internal convention (see makepdf.f90: model%moment=10**(1.5*mw+9.1),
# the N*m/SI form, which only works dimensionally if area/length are
# m/m^2 internally).

_AREA_CONVERSION_TO_KM2 = {
    "km2": 1.0,
    "m2": 1.0e-6,
}


def _convert_area(value: float, from_unit: str, to_unit: str) -> float:
    if from_unit not in _AREA_CONVERSION_TO_KM2:
        raise ValueError(f"Unknown area unit '{from_unit}'. Use 'km2' or 'm2'.")
    if to_unit not in _AREA_CONVERSION_TO_KM2:
        raise ValueError(f"Unknown area unit '{to_unit}'. Use 'km2' or 'm2'.")
    value_km2 = value * _AREA_CONVERSION_TO_KM2[from_unit]
    return value_km2 / _AREA_CONVERSION_TO_KM2[to_unit]


# ---------------------------------------------------------------------------
# Moment <-> moment-magnitude conversion
# ---------------------------------------------------------------------------
#     Mw = (2/3) * (log10(M0) - C)
#
#     where
#         M0 = scalar seismic moment
#         C  = 9.1  if M0 is in N*m   (SI -- matches k223d's own
#                                        makepdf.f90: moment=10**(1.5*mw+9.1))
#         C  = 16.1 if M0 is in dyne*cm (older CGS convention; 1 N*m =
#                                        1e7 dyne*cm, and 16.1-7=9.1, so
#                                        the two forms are consistent)

_MOMENT_CONSTANT = {"N-m": 9.1, "dyne-cm": 16.1}


def mw_to_moment(mw: float, unit: str = "N-m") -> float:
    """Convert moment magnitude Mw to scalar seismic moment M0."""
    if unit not in _MOMENT_CONSTANT:
        raise ValueError(f"Unknown moment unit '{unit}'. Use 'N-m' or 'dyne-cm'.")
    return 10 ** (1.5 * mw + _MOMENT_CONSTANT[unit])


def moment_to_mw(m0: float, unit: str = "N-m") -> float:
    """Convert scalar seismic moment M0 to moment magnitude Mw."""
    if unit not in _MOMENT_CONSTANT:
        raise ValueError(f"Unknown moment unit '{unit}'. Use 'N-m' or 'dyne-cm'.")
    return (2.0 / 3.0) * (math.log10(m0) - _MOMENT_CONSTANT[unit])


# ---------------------------------------------------------------------------
# Coefficient container
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ScalingCoefficients:
    """
    Coefficients for ONE regression direction of ONE category of ONE
    relation. Which direction `a`/`b` apply to (area-on-magnitude vs.
    magnitude-on-area) is determined by which key this object is stored
    under in a relation's `categories` dict -- see ScalingRelation below.

    Attributes
    ----------
    a, b : regression intercept and slope
    area_unit : the area unit THIS regression was fit in -- "km2" or
        "m2". Metadata about the paper, not a user choice; the user's
        requested unit is handled separately in compute().
    magnitude_type : whether `a`/`b` regress against "Mw" directly (most
        modern relations) or against "M0" i.e. log10(seismic moment)
        (some older relations). compute()'s user-facing in/out is
        always Mw regardless -- this only controls the internal
        conversion step.
    moment_unit : only relevant if magnitude_type == "M0" -- "N-m" (SI)
        or "dyne-cm" (older CGS convention).
    sigma : standard deviation of the regression residuals, if reported.
        UNITS DEPEND ON DIRECTION: log10(area) units for an
        area-on-magnitude fit, Mw units for a magnitude-on-area fit.
    m_min, m_max : magnitude range (Mw) the regression was fit over /
        considered valid for (used only to emit a warning, not to block)
    """
    a: float
    b: float
    area_unit: str = "km2"
    magnitude_type: str = "Mw"
    moment_unit: str = "N-m"
    sigma: Optional[float] = None
    m_min: Optional[float] = None
    m_max: Optional[float] = None


# ---------------------------------------------------------------------------
# Base class: every relation implements this interface
# ---------------------------------------------------------------------------

class ScalingRelation:
    """
    Base class for a magnitude <-> area scaling relationship.

    Subclasses set `name`, `citation`, and `categories`: a dict mapping a
    category label specific to THIS paper -> a dict of up to two entries,
    keyed "area_on_magnitude" and/or "magnitude_on_area", each holding
    the ScalingCoefficients actually published for that direction. A
    category needs only the direction(s) the paper actually provides;
    the other direction is synthesized by algebraic inversion (with a
    warning) if and when it's needed.

    Subclasses may override area_from_magnitude / magnitude_from_area
    entirely if the relation isn't a simple log-linear fit (e.g. Hanks &
    Bakun's bilinear form).
    """

    name: str = "unnamed"
    citation: str = ""
    categories: dict[str, dict[str, ScalingCoefficients]] = {}

    _DIRECTIONS = ("area_on_magnitude", "magnitude_on_area")

    @classmethod
    def list_categories(cls) -> list[str]:
        return list(cls.categories.keys())

    @classmethod
    def list_directions(cls, category: str) -> list[str]:
        """Which of the two regression directions this category has a
        published (non-inverted) fit for."""
        return list(cls._get_category(category).keys())

    @classmethod
    def native_area_unit(cls, category: str) -> str:
        entry = cls._get_category(category)
        return next(iter(entry.values())).area_unit

    @classmethod
    def _get_category(cls, category: str) -> dict[str, ScalingCoefficients]:
        if category not in cls.categories:
            raise ValueError(
                f"Relation '{cls.name}' does not support category "
                f"'{category}'. Supported categories: {cls.list_categories()}"
            )
        return cls.categories[category]

    @classmethod
    def _resolve_direction(cls, category: str, want: str):
        """Return (coefficients, was_inverted) for the requested
        direction, falling back to inverting the other direction's fit
        if the wanted one wasn't published."""
        entry = cls._get_category(category)
        if want in entry:
            return entry[want], False
        other = "magnitude_on_area" if want == "area_on_magnitude" else "area_on_magnitude"
        if other in entry:
            return entry[other], True
        raise ValueError(
            f"Relation '{cls.name}', category '{category}' has no "
            "regression coefficients in either direction."
        )

    _DIRECTION_LABELS = {
        "area-on-magnitude": "area from magnitude",
        "magnitude-on-area": "magnitude from area",
    }

    @classmethod
    def _warn_inversion(cls, relation_name: str, category: str, wanted: str, have: str) -> None:
        warnings.warn(
            f"{relation_name} ({category}): no direct fit for "
            f"{cls._DIRECTION_LABELS[wanted]}; using the "
            f"{cls._DIRECTION_LABELS[have]} fit inverted algebraically -- "
            "treat with caution."
        )

    @staticmethod
    def _warn_if_out_of_range(m: float, c: ScalingCoefficients) -> None:
        if c.m_min is not None and c.m_max is not None:
            if not (c.m_min <= m <= c.m_max):
                warnings.warn(
                    f"Magnitude {m:.2f} is outside the fitted range "
                    f"[{c.m_min}, {c.m_max}] for this relation/category; "
                    "result is an extrapolation."
                )

    @classmethod
    def _to_regression_variable(cls, mw: float, c: ScalingCoefficients) -> float:
        """Convert Mw (always the user-facing input) to whatever variable
        this regression was actually fit against."""
        if c.magnitude_type == "Mw":
            return mw
        elif c.magnitude_type == "M0":
            return math.log10(mw_to_moment(mw, unit=c.moment_unit))
        raise ValueError(f"Unknown magnitude_type '{c.magnitude_type}'.")

    @classmethod
    def _from_regression_variable(cls, x: float, c: ScalingCoefficients) -> float:
        """Convert the regression variable back to Mw (always the
        user-facing output)."""
        if c.magnitude_type == "Mw":
            return x
        elif c.magnitude_type == "M0":
            return moment_to_mw(10 ** x, unit=c.moment_unit)
        raise ValueError(f"Unknown magnitude_type '{c.magnitude_type}'.")

    @classmethod
    def area_from_magnitude(cls, magnitude: float, category: str):
        """
        `magnitude` is always Mw. Returns (area_in_native_unit,
        coefficients_used, was_inverted).
        """
        c, inverted = cls._resolve_direction(category, "area_on_magnitude")
        cls._warn_if_out_of_range(magnitude, c)
        x = cls._to_regression_variable(magnitude, c)
        if not inverted:
            # native form: log10(A) = a + b*x
            area = 10 ** (c.a + c.b * x)
        else:
            cls._warn_inversion(cls.name, category, "area-on-magnitude", "magnitude-on-area")
            # only have x = a + b*log10(A) -> invert for log10(A)
            area = 10 ** ((x - c.a) / c.b)
        return area, c, inverted

    @classmethod
    def magnitude_from_area(cls, area: float, category: str):
        """
        `area` must be in this relation's NATIVE unit. Returns
        (Mw, coefficients_used, was_inverted).
        """
        c, inverted = cls._resolve_direction(category, "magnitude_on_area")
        if not inverted:
            # native form: x = a + b*log10(A)
            x = c.a + c.b * math.log10(area)
        else:
            cls._warn_inversion(cls.name, category, "magnitude-on-area", "area-on-magnitude")
            # only have log10(A) = a + b*x -> invert for x
            x = (math.log10(area) - c.a) / c.b
        mw = cls._from_regression_variable(x, c)
        cls._warn_if_out_of_range(mw, c)
        return mw, c, inverted


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_REGISTRY: dict[str, type[ScalingRelation]] = {}


def register_relation(key: str):
    """Class decorator to register a ScalingRelation under a short key."""
    def _wrap(cls: type[ScalingRelation]):
        _REGISTRY[key] = cls
        return cls
    return _wrap


def available_relations() -> dict[str, list[str]]:
    """Return {relation_key: [supported categories]} for discovery."""
    return {key: cls.list_categories() for key, cls in _REGISTRY.items()}


# ---------------------------------------------------------------------------
# Concrete relations
# ---------------------------------------------------------------------------

@register_relation("wells_coppersmith_1994")
class WellsCoppersmith1994(ScalingRelation):
    """
    Wells, D.L. & Coppersmith, K.J. (1994), "New empirical relationships
    among magnitude, rupture length, rupture width, rupture area, and
    surface displacement", BSSA 84(4), 974-1002.

    Subdivided by faulting style. This paper publishes BOTH regression
    directions as separate fits: Table 2A gives area-on-magnitude
    (log10(RA) = a + b*Mw) and magnitude-on-area
    (Mw = a + b*log10(RA)). 
    """
    name = "Wells & Coppersmith (1994)"
    citation = "BSSA 84(4), 974-1002"
    categories = {
        "strike-slip": {
            "area_on_magnitude": ScalingCoefficients(a=-3.42, b=0.90, sigma=0.22, m_min=4.8, m_max=7.9,area_unit="km2"),
            "magnitude_on_area": ScalingCoefficients(a=3.98, b=1.02, sigma=0.23, m_min=4.8, m_max=7.9,area_unit="km2"),
        },
        "reverse": {
            "area_on_magnitude": ScalingCoefficients(a=-3.99, b=0.98, sigma=0.26, m_min=4.8, m_max=7.6,area_unit="km2"),
            "magnitude_on_area": ScalingCoefficients(a=4.33, b=0.90, sigma=0.25, m_min=4.8, m_max=7.6,area_unit="km2"),
        },
        "normal": {
            "area_on_magnitude": ScalingCoefficients(a=-2.87, b=0.82, sigma=0.22, m_min=5.2, m_max=7.3,area_unit="km2"),
            "magnitude_on_area": ScalingCoefficients(a=3.93, b=1.02, sigma=0.25, m_min=5.2, m_max=7.3,area_unit="km2"),
        },
        "all": {
            "area_on_magnitude": ScalingCoefficients(a=-3.49, b=0.91, sigma=0.24, m_min=4.8, m_max=7.9,area_unit="km2"),
            "magnitude_on_area": ScalingCoefficients(a=4.07, b=0.98, sigma=0.24, m_min=4.8, m_max=7.9,area_unit="km2"),
        },
    }


@register_relation("strasser_2010")
class Strasser2010(ScalingRelation):
    """
    Strasser, F.O., Arango, M.C. & Bommer, J.J. (2010), "Scaling of the
    source dimensions of interface and intraslab subduction-zone
    earthquakes with moment magnitude", Seismological Research Letters
    81(6), 941-950.

    Subdivided by subduction-zone rupture position, NOT faulting style.
    Only the area-on-magnitude direction is included here; if this paper
    also published a magnitude-on-area fit and you want it, add a
    "magnitude_on_area" entry to the relevant category below.
    """
    name = "Strasser et al. (2010)"
    citation = "SRL 81(6), 941-950"
    categories = {
        "interface": {"area_on_magnitude": ScalingCoefficients(a=-3.476, b=0.952, sigma=0.304, m_min=6.3, m_max=9.4,area_unit="km2"),
                    "magnitude_on_area": ScalingCoefficients(a=4.441, b=0.846, sigma=0.286, m_min=6.3, m_max=9.4,area_unit="km2")        
        },
        "intraslab": {"area_on_magnitude": ScalingCoefficients(a=-3.225, b=0.890, sigma=0.184, m_min=5.9, m_max=7.8,area_unit="km2"),
                      "magnitude_on_area": ScalingCoefficients(a=4.054, b=0.981, sigma=0.193, m_min=5.9, m_max=7.8,area_unit="km2"),
        },
    }


@register_relation("leonard_2010")
class Leonard2010(ScalingRelation):
    """
    Leonard, M. (2010), "Earthquake fault scaling: self-consistent
    relating of rupture length, width, average displacement and moment
    release", BSSA 100(5A), 1971-1988.

    Taken from Table 6
    """
    name = "Leonard (2010)"
    citation = "BSSA 100(5A), 1971-1988"
    categories = {
        "strike-slip": {"magnitude_on_area": ScalingCoefficients(a=3.99, b=1.0,area_unit="km2")},
        "dip-slip":    {"magnitude_on_area": ScalingCoefficients(a=4.00, b=1.0,area_unit="km2")},
        "interplate": {"magnitude_on_area": ScalingCoefficients(a=4.19, b=1.0,area_unit="km2")},
    }


@register_relation("hanks_bakun_2002")
class HanksBakun2002(ScalingRelation):
    """
    Hanks, T.C. & Bakun, W.H. (2002), "A bilinear source-scaling model for
    M-logA observations of continental earthquakes", BSSA 92(5), 1841-1846.

    Only the "strike-slip" category is registered below (the paper's own
    dataset is predominantly strike-slip); requesting any other category
    raises via the base class's _get_category(). The M-logA relation is
    BILINEAR (the slope changes above a hinge area of 10^2.73 ~ 537 km^2,
    equivalent to Mw 6.71), so this relation overrides the base class
    methods entirely rather than using a simple a + b*x lookup; the
    ScalingCoefficients entry below is kept only so category
    listing/validation works consistently with other relations -- its
    a/b fields aren't used.
    """
    name = "Hanks & Bakun (2002)"
    citation = "BSSA 92(5), 1841-1846"
    categories = {"strike-slip": {"area_on_magnitude": ScalingCoefficients(a=None, b=None, area_unit="km2")}}

    _Mw_HINGE = 6.71  # moment magnitude at the break in slope
    _A_HINGE = 10 ** 2.73  # km^2, area at the break in slope (native unit for this relation) -- self-consistent with _Mw_HINGE via both branches below

    @classmethod
    def magnitude_from_area(cls, area: float, category: str = "strike-slip"):
        c = cls._get_category(category)["magnitude_on_area"]
        if area <= cls._A_HINGE:
            mw = math.log10(area) + 3.98
        else:
            mw = (4.0 / 3.0) * math.log10(area) + 3.07
        return mw, c, False

    @classmethod
    def area_from_magnitude(cls, magnitude: float, category: str = "strike-slip"):
        c = cls._get_category(category)["area_on_magnitude"]
        if magnitude <= cls._Mw_HINGE:
            area = 10 ** (magnitude - 3.98)
        else:
            area = 10 ** ((magnitude - 3.07) * (3.0 / 4.0))
        return area, c, False


@register_relation("_illustrative")
class _Illustrative(ScalingRelation):
    """
    NOT a real published relation -- included only to demonstrate (a) the
    magnitude_type="M0" pathway for relations fit against log10(moment)
    rather than Mw directly, and (b) the automatic direction-selection /
    fallback-with-warning mechanism. Replace or remove once you've added
    real relations that need these features.
    """
    name = "Illustrative relation (placeholder, not a real paper)"
    citation = ""
    categories = {
        # Has both directions -- compute() will use whichever is the
        # statistically correct one for the operation, never inverting.
        "example_both": {
            "area_on_magnitude": ScalingCoefficients(a=-9.5, b=0.66, area_unit="km2", magnitude_type="M0"),
            "magnitude_on_area": ScalingCoefficients(a=14.39, b=1.515, area_unit="km2", magnitude_type="M0"),
        },
        # Has only area-on-magnitude -- asking for magnitude-from-area
        # will trigger the inversion fallback and its warning.
        "example_area_only": {
            "area_on_magnitude": ScalingCoefficients(a=-9.5, b=0.66, area_unit="km2", magnitude_type="M0"),
        },
    }


# ---------------------------------------------------------------------------
# Top-level user-facing dispatch function
# ---------------------------------------------------------------------------

def compute(
    relation: str,
    category: str,
    magnitude: Optional[float] = None,
    area: Optional[float] = None,
    area_unit: str = "m2",
    return_coefficients: bool = False,
):
    """
    Compute area from magnitude, or magnitude from area, using the named
    scaling relation and category. Magnitude is always Mw, in and out.

    Automatically uses whichever regression direction was actually
    published for the operation you're doing (area-on-magnitude when
    solving for area, magnitude-on-area when solving for magnitude);
    only falls back to algebraically inverting the other direction, with
    a warning, if the needed direction wasn't published for that
    category.

    Parameters
    ----------
    relation : registry key, e.g. "wells_coppersmith_1994" (see
        available_relations() for the full list)
    category : the category label THAT RELATION supports, e.g.
        "strike-slip" for Wells & Coppersmith, "interface" for Strasser.
    magnitude : moment magnitude Mw, if solving for area
    area : rupture area, if solving for magnitude -- given in `area_unit`
    area_unit : the unit YOU are working in ("km2" or "m2"), for both the
        `area` argument and the returned area. Defaults to "m2" to match
        k223d's own internal convention.
    return_coefficients : if True, return (value, coefficients) instead
        of just value, where coefficients is the ScalingCoefficients
        object actually used (a, b, sigma, area_unit, etc. -- for
        whichever direction was used) -- useful for logging provenance.

    Exactly one of `magnitude` or `area` must be given.
    """
    if (magnitude is None) == (area is None):
        raise ValueError("Provide exactly one of `magnitude` or `area`.")

    if relation not in _REGISTRY:
        raise ValueError(f"Unknown relation '{relation}'. Available: {list(_REGISTRY)}")
    cls = _REGISTRY[relation]
    native_unit = cls.native_area_unit(category)

    if magnitude is not None:
        area_native, coeffs, _inverted = cls.area_from_magnitude(magnitude, category)
        result = _convert_area(area_native, from_unit=native_unit, to_unit=area_unit)
    else:
        area_native = _convert_area(area, from_unit=area_unit, to_unit=native_unit)
        result, coeffs, _inverted = cls.magnitude_from_area(area_native, category)

    if return_coefficients:
        return result, coeffs
    return result


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print(available_relations())

    # Default output unit is m^2, matching k223d's internal convention
    A = compute("wells_coppersmith_1994", "strike-slip", magnitude=7.0)
    print(f"W&C94 strike-slip, Mw 7.0 -> Area = {A:.3e} m^2 (uses Table 2A fit, no warning)")

    # Round trip using the DIRECT magnitude-on-area fit (Table 2B) -- not
    # the same numbers you'd get by algebraically inverting Table 2A
    M_direct = compute("wells_coppersmith_1994", "strike-slip", area=A, area_unit="m2")
    M_inverted = (math.log10(A * 1e-6) - (-3.42)) / 0.90  # manual inversion of Table 2A, for comparison
    print(f"W&C94 strike-slip, Area -> Mw (Table 2B direct fit) = {M_direct:.4f}")
    print(f"W&C94 strike-slip, Area -> Mw (naive Table 2A inversion) = {M_inverted:.4f}  <- these differ")

    # Strasser only has area-on-magnitude published here -> asking for
    # magnitude from area triggers the fallback + warning
    print("\nExpect a UserWarning below (Strasser has no magnitude-on-area fit):")
    M_strasser = compute("strasser_2010", "interface", area=5000e6, area_unit="m2")
    print(f"Strasser10 interface, Area 5000e6 m^2 -> Mw = {M_strasser:.2f}")

    # Hanks & Bakun bilinear form
    M_hb = compute("hanks_bakun_2002", "strike-slip", area=200, area_unit="km2")
    print(f"\nH&B02, Area 200 km^2 -> Mw = {M_hb:.2f}")

    # Illustrative M0-based relation, both directions available
    A_illust = compute("_illustrative", "example_both", magnitude=7.0)
    M_illust = compute("_illustrative", "example_both", area=A_illust, area_unit="m2")
    print(f"\nIllustrative (both directions, M0-based), Mw 7.0 -> Area = {A_illust:.3e} m^2 -> Mw = {M_illust:.3f}")

    # k223d's own Mw<->moment formula (makepdf.f90 line 208) reproduced exactly
    mw_test = 7.0
    m0_test = mw_to_moment(mw_test, unit="N-m")
    print(f"\nMw {mw_test} -> M0 = {m0_test:.3e} N*m (matches k223d's 10**(1.5*mw+9.1))")

    # Getting the coefficients back alongside the result
    A3, coeffs = compute(
        "wells_coppersmith_1994", "strike-slip", magnitude=7.0, return_coefficients=True
    )
    print(f"\nArea = {A3:.3e} m^2, using a={coeffs.a}, b={coeffs.b}, sigma={coeffs.sigma}")
