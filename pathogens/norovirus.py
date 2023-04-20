from pathogen_properties import *

background = """Norivirus in a GI infection, mostly spread through personal
contact."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=142786,
)


variables = {
    "cases": IncidenceAbsolute(
        annual_infections=20e6,
        confidence_interval=(19e6, 21e6),
        country="United States",
        source="https://www.cdc.gov/norovirus/trends-outbreaks/burden-US.html",
    ),
    "us_population": Population(
        people=333_287_557,
        country="United States",
        source="https://www.census.gov/quickfacts/fact/table/US/PST045221",
        start_date="2022-07-01",
        end_date="2022-07-01",
    ),
    "shedding_duration": SheddingDuration(
        days=2,
        confidence_interval=(1, 3),
        source="https://www.mayoclinic.org/diseases-conditions/norovirus/symptoms-causes/syc-20355296",
    ),
    "rothman_period_outbreaks": Number(
        number=13 / 5,  # monthly outbreaks
        country="United States",
        source="https://www.cdc.gov/norovirus/reporting/norostat/data-table.html",
    ),
    "normal_year_outbreaks": Number(
        number=1246 / 12,  # monthly outbreaks
        country="United States",
        source="https://www.cdc.gov/norovirus/reporting/norostat/data-table.html",
    ),
}


def estimate_prevalences():
    return {
        "rothman_2021_period": variables["cases"]
        .to_rate(variables["us_population"])
        .to_prevalence(variables["shedding_duration"])
        .scale(
            variables["rothman_period_outbreaks"].per(
                variables["normal_year_outbreaks"]
            )
        )
    }
