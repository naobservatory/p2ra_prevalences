import dataclasses

from pathogen_properties import *

background = """Adeno-associated virus 2 is a small replication-defective
virus"""

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(10804),
)


def estimate_prevalences() -> list[Prevalence]:
    stub_estimate = Prevalence(
        # Completely made up, checking in stub for model testing.
        infections_per_100k=50_000,
        active=Active.LATENT,
        country="United States",
    )
    return [
        # This pathogen should be close to constant, so extrapolate to 2020 and
        # 2021.
        dataclasses.replace(stub_estimate, date_source=Variable(date="2020")),
        dataclasses.replace(stub_estimate, date_source=Variable(date="2021")),
    ]


def estimate_incidences():
    return []
