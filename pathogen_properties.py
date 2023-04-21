from enum import Enum
from dataclasses import dataclass
from typing import Optional

# Enums, short for enumerations, are a data type in Python used to represent a set of named values,
# which are typically used to define a set of related constants with unique names.

# TODO: Potentially create enum for country, state, etc. to ensure consistency in naming


class NAType(Enum):
    DNA = "DNA"
    RNA = "RNA"


class Enveloped(Enum):
    ENVELOPED = "enveloped"
    NON_ENVELOPED = "non_enveloped"


class VariableType(Enum):
    MEASUREMENT = "measurement"
    EXTERNAL_ESTIMATE = "external_estimate"
    NAO_ESTIMATE = "nao_estimate"


@dataclass(kw_only=True)
class PathogenChars:
    na_type: NAType
    enveloped: Enveloped
    taxid: int


@dataclass(kw_only=True)
class Variable:
    """An external piece of data"""

    source: Optional[str] = None
    country: Optional[str] = None
    state: Optional[str] = None
    county: Optional[str] = None
    number_of_participants: Optional[int] = None
    confidence_interval: Optional[tuple[float, float]] = None
    methods: Optional[str] = None
    # Either supply date, or start_date and end_date.
    # Dates can be any of: YYYY, YYYY-MM, or YYYY-MM-DD.
    date: Optional[str] = None
    start_date: Optional[str] = None
    end_date: Optional[str] = None

    # Remember to recursively consider each input's inputs if defined.
    inputs: Optional[list["Variable"]] = None


@dataclass(kw_only=True)
class Population(Variable):
    """A number of people"""

    people: float


@dataclass(kw_only=True)
class Scalar(Variable):
    scalar: float


@dataclass(kw_only=True)
class Prevalence(Variable):
    """What fraction of people have this pathogen at some moment"""

    infections_per_100k: float

    def scale(self, scalar: Scalar) -> "Prevalence":
        return Prevalence(
            infections_per_100k=self.infections_per_100k * scalar.scalar,
            inputs=[self, scalar],
        )


@dataclass(kw_only=True)
class PrevalenceAbsolute(Variable):
    """How many people had this pathogen at some moment"""

    infections: float

    def to_rate(self, population: Population) -> Prevalence:
        return Prevalence(
            infections_per_100k=self.infections * 100000 / population.people,
            inputs=[self, population],
        )


@dataclass(kw_only=True)
class SheddingDuration(Variable):
    days: float


@dataclass(kw_only=True)
class Number(Variable):
    """Generic number.  Use this for weird one-off things

    If some concept is being used by more than two estimates it should get a
    more specific Variable subclass."""

    number: float

    def per(self, other: "Number") -> Scalar:
        return Scalar(scalar=self.number / other.number, inputs=[self, other])


@dataclass(kw_only=True)
class IncidenceRate(Variable):
    """What fraction of people get this pathogen annually"""

    annual_infections_per_100k: float

    def to_prevalence(self, shedding_duration: SheddingDuration) -> Prevalence:
        return Prevalence(
            infections_per_100k=self.annual_infections_per_100k
            * shedding_duration.days
            / 365,
            inputs=[self, shedding_duration],
        )


@dataclass(kw_only=True)
class IncidenceAbsolute(Variable):
    """How many people get this pathogen annually"""

    annual_infections: float

    def to_rate(self, population: Population) -> IncidenceRate:
        return IncidenceRate(
            annual_infections_per_100k=self.annual_infections
            * 100000
            / population.people,
            inputs=[self, population],
        )
